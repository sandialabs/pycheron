#####################################################################################
# Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
# certain rights in this software.
#####################################################################################
# NOTICE:
# For five (5) years from 10/21/2019 the United States Government is granted for
# itself and others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and perform publicly and
# display publicly, by or on behalf of the Government. There is provision for the
# possible extension of the term of this license. Subsequent to that period or any
# extension granted, the United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this
# data to reproduce, prepare derivative works, distribute copies to the public,
# perform publicly and display publicly, and to permit others to do so. The specific
# term of the license can be identified by inquiry made to National Technology and
# Engineering Solutions of Sandia, LLC or DOE. NEITHER THE UNITED STATES GOVERNMENT,
# NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR NATIONAL TECHNOLOGY AND ENGINEERING
# SOLUTIONS OF SANDIA, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR
# IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
# USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. Any licensee of this software
# has the obligation and responsibility to abide by the applicable export control laws,
# regulations, and general prohibitions relating to the export of technical data.
# Failure to obtain an export control license or other authority from the Government
# may result in criminal liability under U.S. laws.
# (End of Notice)
####################################################################################

__all__ = ["transferFunctionMetric"]

import scipy.signal as sig
import numpy as np
from math import floor, log, pi
from pycheron.psd.crossSpectrum import crossSpectrum
from pycheron.psd.McNamaraBins import McNamaraBins
from obspy.clients.iris import Client
from pycheron.util.logger import Logger


def transferFunctionMetric(tr1, tr2, logger=None, database=None):
    """
    Calculates metrics that assess the relationship between two SNCLs with the same network, station and channel but
    separate locations.

    :param tr1: (trace) obspy trace object 1
    :type tr1: obspy.core.trace.Trace
    :param tr2: (trace) obspy trace object 2
    :type tr2: obspy.core.trace.Trace
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database

    :return: (dict) - dictionary with the following keys and values:

                    * snclq (`str`)
                    * start_time (`str`)
                    * end_time (`str`)
                    * metric_name (`str`)
                    * gain_ratio (`float`)
                    * phase_diff (`float`)
                    * ms_coherence (`float`)

    :rtype: dict

    When seismometers are working properly, the transfer function amplitude and phase will match similar values
    calculated from the instrument responses.

    This function calculates the transfer function from data in the incoming streams (traces really). Response
    information is then obtained from the evalresp web service [#]_

    .. note::

        Seismic streams passed to transferFunctionMetric() must have the same network, station, and channel and must
        cover the same time range. The two channels should also have values of azimuth and dip within 5 degrees of each
        other. If sampling rates differ, and one is a multiple of the other, the stream with the higher sampling rate will
        be decimated to match the lower sampling frequency. The output dictionary generated for these two-channel metrics
        will have a SNCL code of the form: Net.Sta.Loc1:Loc2.Chan.

    **Algorithm Steps**

    1) Calculate complex transfer function:
         * Request 1 hour of data (noise is ideal, so any consistent time is fine; only need this monthly or weekly)
         * If channel sample rates differ, decimate to make them identical
         * Detrend
         * Demean
         * Compute complex cross-spectrum of traces x and y --> Pxx, Pxy, Pyy (same procedure as PSD)
         * Retrieve stored PSD Pxx(f)
    2) Calculate transfer function values:
         * Txy(f) = Pxy(f) / Pxx(f)
         * dataGain <- Mod(Txy)
         * dataPhase <- Arg(Txy)
         * Save avgDataGain (amplitude) and avgDataphase values for T=5-7s
    3) Calculate the corresponding response amplitude ratio and phase difference for x and y
         * Request responses for x and y for that 1 hour
         * Calculate respGain = respGainy(f) / respGainx(f)
         * Calculate respPhase = respPhasey(f) - respPhasex(f)
         * Save avgRespGain and avgRespPhase values for T = 5-7s
    4) Calculate metrics:
         * Calculate (data/resp) gain ratio for two saved gain averages:
           gain_ratio = avgDataGain/avgRespGain (reasonableness of cross-spectral amplitudes between st1 and
           st2)
         * Calculate (data-resp) phase difference from the two saved phase averages:
           phase_diff = avgDataPhase - avgRespPhase (reasonableness of cross-spectral phases between st1 and
           st2)
         * Calculate magnitude squared coherence of x and y:

           * Retrieve stored PSD for Pxx and Pyy
           * ms_coherence = abs(Pxy)^2 / (Pxx*Pyy) (mean square coherence between st1 and st2)
           * Save the average MSCoherence for T=5-7s

    These values can be interpreted as follows:

    * Whenever ms_coherence ~= 1.0, properly functioning seismometers should have:

      * gain ratio ~= 1.0
      * phase_diff < 10.0 (degrees)

    **Example**

    .. code-block:: python

        # Initialize IRIS client
        client = Client("IRIS")

        #Define start/end times
        starttime = UTCDateTime("2011-05-01T01:00:00.000")
        endtime = starttime + 3600

        # Grab data from IRIS web client server with specified starttime/endtime and SNCL
        st1 = client.get_waveforms("CI","PASC","00","BHZ", starttime,endtime)
        st2 = client.get_waveforms("CI","PASC","10","BHZ", starttime, endtime)

        # Grab out traces
        tr1 = st1[0]
        tr2 = st2[0]

        # Calculate transfer function then print
        TF = transferFunctionMetric(tr1,tr2)
        print TF
        >>> {'metric_name': 'transferFunction', 'start_time': "2011-05-01T01:00:00.019500Z", 'snclq': 'PASC:CI:BHZ:00:BHZ:10', 'end_time': "2011-05-01T01:59:59.994500Z", 'ms_coherence': 0.9999875790954715, 'gain_ratio': 0.9058981533571073, 'phase_diff': 0.13553317118827302}

    .. rubric:: References

    .. [#] http://service.iris.edu/irisws/evalresp/1/

    """

    # --------------Sanity Checks and data re-sampling-----------------------------------------------
    if logger == None:
        logger = Logger(None)

    # First though let's define default span for crossSpectrum calculation
    spans = np.array([3, 5])

    # Get start and end times for each trace
    start1 = tr1.stats.starttime
    end1 = tr1.stats.endtime

    start2 = tr2.stats.starttime
    end2 = tr2.stats.endtime

    # Checking temporal extent of traces match, if not exit
    if start2 - start1 > 1 or start2 - start1 < 0:
        logger.warn(
            "transferFunctionMetric(): Start times do not match for "
            "traces %s:%s" % (str(tr1.get_id()), str(tr2.get_id()))
        )
        return

    if end2 - end1 > 1 or end2 - end1 < 0:
        logger.warn(
            "transferFunctionMetric(): End times do not match for "
            "traces {}:{}".format(str(tr1.get_id()), str(tr2.get_id()))
        )
        return

    # Checking sampling rates > 1
    if tr1.stats.sampling_rate < 1:
        logger.warn(
            "transferFunctionMetric(): {} has a sampling rate < 1".format(
                str(tr1.get_id())
            )
        )
        return

    if tr2.stats.sampling_rate < 1:
        logger.warn(
            "transferFunctionMetric(): {} has a sampling rate < 1".format(
                (str(tr2.get_id()))
            )
        )
        return

    # Check if sampling rates differ, make sampling rate minimum of both traces if differ
    sr1 = round(tr1.stats.sampling_rate)
    sr2 = round(tr2.stats.sampling_rate)
    sampling_rate = min(sr1, sr2)

    # Check sample rate for trace 1 is multiple of trace 2. If different decimate so they are the same. If they are the
    # same assign data1 = tr1.data (e.g., decimation not necessary, so data is just original data)
    if sr1 > sampling_rate:
        if sr1 % sampling_rate != 0:
            logger.warn(
                "transferFunctionMetric(): Sampling Rates are not multiples of "
                "each other {s1}:{s2}".format(
                    s1=str(tr1.get_id()), s2=str(tr2.get_id())
                )
            )
            return
        increment = round(sr1 / sampling_rate)
        data1 = sig.decimate(tr1.data, int(increment))
    else:
        data1 = tr1.data

    # Check sample rate for trace 2 is multiple of trace 1. If different decimate so they are the same. If they are the
    # same assign data2 = tr2.data (e.g., decimation not necessary, so data is just original data)
    if sr2 > sampling_rate:
        if sr2 % sampling_rate != 0:
            logger.warn(
                "transferFunctionMetric(): Sampling Rates are not multiples of "
                "each other %s:%s" % (str(tr1.get_id()), str(tr2.get_id()))
            )
            return

        increment = round(sr2 / sampling_rate)
        data2 = sig.decimate(tr2.data, int(increment))
    else:
        data2 = tr2.data

    # Check for valid data everywhere after potential resampling
    if np.isnan(data1).any() or np.isnan(data2).any():
        logger.warn(
            "transferFunctionMetric(): NaN values generated during "
            "re-sampling %s:%s" % (str(tr1.get_id()), str(tr2.get_id()))
        )
        return

    # ----------------Spectral Analysis-----------------------------------------

    # Choose McNamara frequencies (lo,hi) based on channel band code (e.g, sampling rate). For more information see:
    # http://www.iris.edu/manuals/SEED_appA.htm

    channel = tr1.stats.channel
    if channel.startswith("L"):
        loFreq = 0.001
        hiFreq = round(0.5 * tr1.stats.sampling_rate)
    elif channel.startswith("M"):
        loFreq = 0.0025
        hiFreq = round(0.5 * tr1.stats.sampling_rate)
    else:
        loFreq = 0.005
        hiFreq = round(0.5 * tr1.stats.sampling_rate)

    # Set an alignment frequency from which octaves will be generated
    alignFreq = 0.1

    # Truncate segment to nearest power of 2 samples
    pow2 = floor(log(len(data1), 2))
    truncateLen = int(2 ** pow2)

    # initializing summation arrays for truncated trace data divided into 13 segments
    specSum1 = np.zeros(1)
    specSum2 = np.zeros(1)
    pxxSum = np.zeros(1)
    pyySum = np.zeros(1)
    pxySum = np.zeros(1)
    cohSum = np.zeros(1)
    phaSum = np.zeros(1)

    # Divide the truncated trace data into 13 segments with 75% overlap. Calculate crossSpectrum for each segment.
    for i in range(0, 13):
        first = int(i * (truncateLen / 16))
        last = int((i + 4) * truncateLen / 16)
        d1 = data1[first:last]
        d2 = data2[first:last]
        # Calculate cross spectrum for each segment
        freq, spec1, spec2, coh, phase, pxx, pyy, pxy = crossSpectrum(
            d1,
            d2,
            sampling_rate,
            spans=spans,
            taper=0.1,
            pad=0,
            demean=False,
            detrend=True,
            logger=logger,
        )
        # csd and periodogram returns 'two-sided' spectra and needs to be multiplied
        # by 2.0 when the time series being evaluated is Real rather than Complex.
        # Please see these two excellent sources, especially the 'psd' vignette:
        #   http://www.stanford.edu/class/ee262/software/signal_tb.pdf
        #   http://cran.r-project.org/web/packages/psd/vignettes/normalization.pdf

        # summing metrics - sum each segment to obtain final spec1, spec2, pxx, pyy, pxy, coherence, and phase
        specSum1 = specSum1 + 2 * spec1
        specSum2 = specSum2 + 2 * spec2
        pxxSum = pxxSum + pxx
        pyySum = pyySum + pyy
        pxySum = pxySum + pxy
        cohSum = cohSum + coh
        phaSum = phaSum + phase

    # Get average of the sums for each of above variables

    pxx = pxxSum / 13
    pxy = pxySum / 13
    coh = cohSum / 13
    phase = phaSum / 13

    # Get McNamara binned version of each variable
    # Note, when run on components with imaginary components throws warning: ComplexWarning: Casting complex
    # values to real discards the imaginary part dfp[i] = (sum(pxx[idx_lo:idx_hi]) / (idx_hi - idx_lo)). Pxx, Pyy
    # imaginary components were small enough (*e-10/e-11) that their imaginary component really didn't affect answers,
    # though imaginary component of pxy affects dataAmp answer. In order to counter this difference, computed real and
    # imaginary components separately then added them back together afterward. Final output result agrees with R code
    # output.

    f, PXXR = McNamaraBins(freq, np.real(pxx), loFreq, hiFreq, alignFreq)
    f, PXYR = McNamaraBins(freq, np.real(pxy), loFreq, hiFreq, alignFreq)
    f, coh = McNamaraBins(freq, coh, loFreq, hiFreq, alignFreq)
    f, phase = McNamaraBins(freq, phase, loFreq, hiFreq, alignFreq)
    f, PXXI = McNamaraBins(freq, np.imag(pxx), loFreq, hiFreq, alignFreq)
    f, PXYI = McNamaraBins(freq, np.imag(pxy), loFreq, hiFreq, alignFreq)
    pxx = PXXR + PXXI * 1j
    pxy = PXYR + PXYI * 1j

    # ----------------Transfer Function Metrics----------------------------------------

    # Transfer function average amplitude and phase over periods of 5-7 seconds.
    # Need to convert phase from radians to positive degrees for comparison with evalresp.
    dataAmp = abs(pxy / pxx)
    dataPhase = phase * 180 / pi

    # Calculate average values for periods in the 5-7 second range for amp, phase, and coherence.
    # First find indices that fit within period range.
    period = 1 / f
    indices = []
    for i in range(len(period)):
        if 5 <= period[i] <= 7:
            indices.append(i)

    avgDataAmp = np.mean(dataAmp[indices])
    avgDataPhase = np.mean(dataPhase[indices])
    avgCoherence = np.mean(coh[indices])

    # Grab out network, station, channel,location information from traces for evalresp
    network1 = tr1.stats.network
    network2 = tr2.stats.network

    station1 = tr1.stats.station
    station2 = tr2.stats.station

    channel1 = tr1.stats.channel
    channel2 = tr2.stats.channel

    location1 = tr1.stats.location
    location2 = tr2.stats.location

    # Get min/max, len(f) and define units for evalresp
    minfreq = min(f)
    maxfreq = max(f)
    nfreq = len(f)
    units = "def"

    # calling iris client server to obtain evalresp fap response files
    client = Client(timeout=30)
    iris1 = client.evalresp(
        network=network1,
        station=station1,
        location=location1,
        channel=channel1,
        time=start1,
        minfreq=minfreq,
        maxfreq=maxfreq,
        nfreq=nfreq,
        units=units,
        output="fap",
    )

    iris2 = client.evalresp(
        network=network2,
        station=station2,
        location=location2,
        channel=channel2,
        time=start2,
        minfreq=minfreq,
        maxfreq=maxfreq,
        nfreq=nfreq,
        units=units,
        output="fap",
    )

    # Calculate the corresponding response amplitude ratio and phase difference for x and y. In order to compare against
    # pxy/pxx need iris2(amp)/iris1(amp)
    respAmp = iris2[:, 1] / iris1[:, 1]
    respPhase = iris1[:, 2] - iris2[:, 2]

    # Get indices that are within 5-7 period range for evalresp responses
    period = 1 / iris1[:, 0]
    index = []
    for i in range(len(period)):
        if 5 <= period[i] <= 7:
            index.append(i)

    # Calculate average values for periods in the 5-7 second range for evalresp responses (amp, phase)
    avgRespAmp = np.mean(respAmp[index])
    avgRespPhase = np.mean(respPhase[index])

    # Calculate ratio of transfer function average amplitude/phase to response average amplitude/phase over 5-7s periods
    dataRespGainRatio = avgDataAmp / avgRespAmp
    dataRespPhaseDiff = avgDataPhase - avgRespPhase

    # Calculate dataRespPhaseDiff
    if abs(dataRespPhaseDiff) > 180:
        dataRespPhaseDiffMag = 360 - dataRespPhaseDiff
    else:
        dataRespPhaseDiffMag = dataRespPhaseDiff

    # Create combo snclq for both traces, then return output
    snclq = "{0}:{1}".format(tr1.get_id(), tr2.get_id())

    d = {
        "snclq": snclq,
        "start_time": start1.isoformat(),
        "end_time": end1.isoformat(),
        "metric_name": "transferFunctionMetric",
        "gain_ratio": dataRespGainRatio,
        "phase_diff": dataRespPhaseDiffMag,
        "ms_coherence": avgCoherence,
    }

    if database is not None:
        database.insert_metric(d)

    return d
