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

__all__ = ["crossCorrMetric"]

import numpy as np
from pycheron.psd.noise.deadChannel import DDT, isDC
from pycheron.util.logger import Logger
from obspy.signal import cross_correlation, filter


def crossCorrMetric(
    tr1,
    tr2,
    maxLagSecs=10,
    filt="lowpass",
    freqmin=1.0,
    freqmax=100.0,
    corners=4,
    zerophase=False,
    maxorder=None,
    ba=False,
    freq_passband=False,
    logger=None,
    database=None,
):
    """
    Cross-correlates the current data with the neighboring data.

    :param tr1: obspy trace object 1
    :type tr1: `obspy.core.trace.Trace`
    :param tr2: obspy trace object 2
    :type tr2: `obspy.core.trace.Trace`
    :param maxLagSecs: Maximum number of seconds of lag to use
    :type maxLagSecs: int
    :param filt: Choose from the following filter options: `"bandpass"`, `"bandstop"`, `"envelope"`, `"highpass"`, `"integer_decimation"`, `"lowpass"`, `"lowpass_cheby_2"`. See Filter_Options_ .
    :type filt: str
    :param freqmin: Depending on the filter utilized, either pass band low corner frequency (bandpass, bandstop), filter corner frequency (lowpass), or the frequency above which signals are attenuated with 95 dB
    :type freqmin: int
    :param freqmax: Depending on the filter utilized, either passband high corner frequency (bandpass, bandstop), or filter corner frequency (highpass)
    :type freqmax: int
    :param corners:  Filter corners/order
    :type corners: int
    :param zerophase: If True, applies filter once forward and once backwards. This results in twice the filter order but zero phase shift in the resulting filtered trace.
    :type zerophase: bool
    :param decimation_factor: Integer decimation factor
    :type decimation_factor: int
    :param maxorder: Maximal order of designed cheby2 filter (option for lowpass_cheby2 only)
    :type maxorder: int
    :param ba: If True, returns only the filter coefficients (b,a) instead of filtering. (option for lowpass_cheby2 only)
    :type ba: bool
    :param freq_passband: If True, return additionally to the filtered data, then iteratively determine the pass band frequency. (option for lowpass_cheby2 only)
    :type freq_passband: bool
    :param logger: If using a logger, (you must create one using the util.logger class)
    :type logger: pycheron.util.logger.Logger
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database

    :return: Dictionary containing the following keys and types:

             * snclq (`str`)
             * snclq1_starttime (`str`)
             * snclq1_endtime (`str`)
             * snclq2_starttime (`str`)
             * snclq2_endtime (`str`)
             * peak_correlation (`float`)
             * peak_lag (`float`)
             * metric_name (`str`)

    :rtype: dict

    Finds the min and max of the cross-correlation function of these, saves the value with the largest absolute value (preserve the sign). Preserves the lag in seconds to the min or max value saved above.
    Finds the difference in traveltime between the 2 traces: (current - neighbor). Adjusts the lag by subtracting this difference (preserving sign) from the lag. In short, calculates the maximum absolute correlation (polarity_check
    (peak correlation)) and lg at maximum correlation (timing_drift, peak lag) associated with two traces.

    .. _Filter_Options:

    **Filter Options**

    * "bandpass" - Butterworth-Bandpass Filter
            *Parameters Required*

            * **freqmin** (`int`) - Pass band low corner frequency
            * **freqmax** (`int`) - Pass band high corner frequency
            * **df** (`int`) - sampling rate in Hz (get from tr.stats.sampling_rate so not needed as input param)
            * **corners** (`int`) - Filter corners/order.
            * **zerophase** (`bool`) - If True, applies filter once forward and once backwards. This results in twice the filter order but zero phase shift in the resulting filtered trace.

    * "bandstop" - Butterworth-Bandstop Filter
            *Parameters Required*

            * **freqmin** (`int`) - Pass band low corner frequency
            * **freqmax** (`int`) - Pass band high corner frequency
            * **df** (`int`) - sampling rate in Hz (get from tr.stats.sampling_rate so not needed as input param)
            * **corners** (`int`) - Filter corners/order.
            * **zerophase** (`bool`) - If True, applies filter once forward and once backwards. This results in twice the filter order but zero phase shift in the resulting filtered trace.

    * "envelope" - Envelope function
        *Parameters Required*

        * **None**

    *  "highpass" - Butterworth-Highpass Filter
            *Parameters Required*

            * **freqmax** (`int`) - Pass band high corner frequency
            * **df** (`int`) - sampling rate in Hz (get from tr.stats.sampling_rate so not needed as input param)
            * **corners** (`int`) - Filter corners/order.
            * **zerophase** (`bool`) - If True, applies filter once forward and once backwards. This results in twice the filter order but zero phase shift in the resulting filtered trace.

    *   "integer_decimation" - Downsampling by applying a simple integer decimation
            *Parameters Required*

            * **decimation_factor** (`int`) - Integer decimation factor

    *    "lowpass" - Butterwoth-Lowpass Filter
            *Parameters Required*

            * **freqmin** (`int`) - Pass band low corner frequency
            * **df** (`int`) - sampling rate in Hz (get from tr.stats.sampling_rate so not needed as input param)
            * **corners** (`int`) - Filter corners/order.
            * **zerophase** (`bool`) - If True, applies filter once forward and once backwards. This results in twice the filter order but zero phase shift in the resulting filtered trace.

    *    "lowpass_cheby_2" - Cheby2- Lowpass Filter
            *Parameters Required*

            * **freqmin** (`int`) - Pass band low corner frequency
            * **maxorder** (`int`) - Maximal order of designed cheby2 filter
            * **ba** ('bool') - If True, returns only the filter coefficients (b,a) instead of filtering. (DEFAULT = False)
            * **freq_passband** (`bool`) - If True, return additionally to the filtered data, then iteratively determine the pass band frequency.

    **Examples**

    .. code-block:: python

        # Initialize IRIS client
        client = Client("IRIS")

        # Define start/end time and then get same signal shifted by 3 seconds
        starttime = UTCDateTime("2013-11-12T07:09:45.000")
        endtime = starttime + 600
        st1 = client.get_waveforms("NM","SLM","00","BHZ",starttime, endtime)
        st2 = client.get_waveforms("NM","SLM","00","BHZ",starttime+3, endtime+3)

        #Get trace data for each stream
        tr1 = st1[0]
        tr2 = st2[0]

        #Cross-correlate
        CC = crossCorrMetric(tr1,tr2,freqmin=0.1, corners = 2)
        print(CC)
        >>> {'peak_correlation': 0.3450079442293356, 'snclq1_starttime': UTCDateTime(2013, 11, 12, 7, 9, 45), 'snclq2_endtime': UTCDateTime(2013, 11, 12, 7, 19, 48), 'snclq2': u'NM.SLM.00.BHZ.M', 'snclq1': u'NM.SLM.00.BHZ.M', 'peak_lag (s)': 0.008580576472234229, 'snclq2_starttime': UTCDateTime(2013, 11, 12, 7, 9, 48), 'metric_name': 'crossCorrMetric', 'snclq1_endtime': UTCDateTime(2013, 11, 12, 7, 19, 45)}
    """

    if logger == None:
        logger = Logger(None)
    # -----------Compatibility Checks and Data Resampling------------------------------------------

    # Demean and detrend both traces (no taper), but first make a copy of the trace so it isn't permanently changed
    trC1 = tr1.copy()
    trC1 = DDT(trC1, True, True, 0)
    trC2 = tr2.copy()
    trC2 = DDT(trC2, True, True, 0)

    # Deal with potentially different sampling rates
    if trC1.stats.sampling_rate < 1 or trC2.stats.sampling_rate < 1:
        sr1 = trC1.stats.sampling_rate
        sr2 = trC2.stats.sampling_rate
    else:
        sr1 = round(float(trC1.stats.sampling_rate))
        sr2 = round(float(trC2.stats.sampling_rate))

    # Get min of both sampling rates
    sampling_rate = min(sr1, sr2)

    # Check if sampling rates are not multiples of each other, first for trace1, then trace2
    if sr1 > sampling_rate:
        if sr1 % sampling_rate > 0:
            logger.error(
                "crossCorrMetric(): sampling rates are not multiples of each other "
                "%s = %s, %s = %s"
                % (str(trC1.get_id()), str(sr1), str(trC2.get_id()), str(sr2))
            )
            return
        increment = round(sr1 / sampling_rate)
        d1 = trC1.decimate(int(increment))
    else:
        d1 = trC1.data

    if sr2 > sampling_rate:
        if sr2 % sampling_rate > 0:
            logger.error(
                "crossCorrMetric(): sampling rates are not multiples of each other "
                "%s = %s, %s = %s"
                % (str(trC1.get_id()), str(sr1), str(trC2.get_id()), str(sr2))
            )
            return
        increment = round(sr2 / sampling_rate)
        d2 = trC2.decimate(int(increment))
    else:
        d2 = trC2.data

    # Check that we have valid data everywhere
    if np.isnan(d1.any()) is True or np.isnan(d2.any()) is True:
        logger.error("crossCorrMetric(): NaN values generated during resampling")
        return

    # Check for flatlined data
    if isDC(trC1):
        logger.error(
            "crossCorrMetric(): %s has one unique sample value (flatlined). Standard deviation is zero,"
            "cross-correlation is undefined" % (str(trC1.get_id()))
        )
        return

    if isDC(trC2):
        logger.error(
            "crossCorrMetric(): %s has one unique sample value (flatlined). Standard deviation is zero,"
            "cross-correlation is undefined" % (str(trC2.get_id()))
        )
        return

    # -----------Correlation Function------------------------------------------

    # Applying low pass or other filter to both traces
    if filt == "bandpass":
        d1 = filter.bandpass(
            d1, freqmin, freqmax, trC1.stats.sampling_rate, corners, zerophase
        )
        d2 = filter.bandpass(
            d2, freqmin, freqmax, trC2.stats.sampling_rate, corners, zerophase
        )
    elif filt == "bandstop":
        d1 = filter.bandstop(
            d1, freqmin, freqmax, trC1.stats.sampling_rate, corners, zerophase
        )
        d2 = filter.bandstop(
            d2, freqmin, freqmax, trC2.stats.sampling_rate, corners, zerophase
        )
    elif filt == "envelope":
        d1 = filter.envelope(d1)
        d2 = filter.envelope(d2)
    elif filt == "highpass":
        d1 = filter.highpass(d1, freqmax, trC1.stats.sampling_rate, corners, zerophase)
        d2 = filter.highpass(d2, freqmax, trC2.stats.sampling_rate, corners, zerophase)
    elif filt == "integer_decimation":
        d1 = filter.envelope(d1)
        d2 = filter.envelope(d2)
    elif filt == "lowpass":
        d1 = filter.lowpass(d1, freqmin, trC1.stats.sampling_rate, corners, zerophase)
        d2 = filter.lowpass(d2, freqmin, trC2.stats.sampling_rate, corners, zerophase)
    elif filt == "lowpass_cheby_2":
        d1 = filter.lowpass_cheby_2(
            d1, freqmin, trC1.stats.sampling_rate, maxorder, ba, freq_passband
        )
        d2 = filter.lowpass_cheby_2(
            d2, freqmin, trC2.stats.sampling_rate, maxorder, ba, freq_passband
        )

    # Calculate cross-correlation
    lagMax = int(sampling_rate * maxLagSecs)
    xcorr = cross_correlation.correlate(d1, d2, lagMax)

    # Find strongest correlation
    peak_correlation = xcorr[1]

    # Get lag associated with the strongest correlation (already provided by cross_correlation.xcorr), then get peak lag
    lagPoints = xcorr[0]
    peak_lag = lagPoints * 1 / sampling_rate

    # Create output metric information
    # if tr is not from mseed, it will not have dq flag
    try:
        snclq1 = tr1.get_id() + "." + tr1.stats.mseed.dataquality
        snclq2 = tr2.get_id() + "." + tr2.stats.mseed.dataquality
    except AttributeError:
        snclq1 = tr1.get_id()
        snclq2 = tr2.get_id()

    d = {
        "snclq": snclq1 + ":" + snclq2,
        "snclq1_starttime": tr1.stats.starttime.isoformat(),
        "snclq1_endtime": tr1.stats.endtime.isoformat(),
        "snclq2_starttime": tr2.stats.starttime.isoformat(),
        "snclq2_endtime": tr2.stats.endtime.isoformat(),
        "peak_correlation": peak_correlation,
        "peak_lag": peak_lag,
        "metric_name": "crossCorrMetric",
    }

    if database is not None:
        database.insert_metric(d)

    return d
