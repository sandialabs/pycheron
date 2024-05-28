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

__all__ = ["psdMetric"]

import numpy as np
from pycheron.psd.psdList import psdList
from pycheron.psd.psdStatistics import psdStatistics
from pycheron.sigpro.unHistogram import unHistogram
import pandas as pd
from collections import OrderedDict
from pycheron.metrics.basicStatsMetric import basicStatsMetric
import warnings
from pycheron.util.logger import Logger
from pycheron.util.masks import samples2time
from pycheron.db.sqllite_db import Database

warnings.filterwarnings("ignore")


def psdMetric(
    st,
    expLoPeriod=None,
    expHiPeriod=100,
    linLoPeriod=None,
    linHiPeriod=100,
    evalresp=None,
    generateMasks=False,
    dcExpThreshold=0.3,
    dcExpThresholdHour=0.25,
    pctBelowNoiseThreshold=20,
    pctAboveNoiseThreshold=20,
    rmsThreshold=50000,
    dcLinThreshold=2,
    dcLinThresholdHour=2,
    pctBelowNoiseThresholdRESP=90,
    pctAboveNoiseThresholdRESP=90,
    processes=6,
    logger=None,
    masksByTime=True,
    byHourOn=False,
    database_config=None,
    session=None,
    overwrite=True,
):
    """
    Performs a power spectral density analysis (based on McNamara method) on seismic traces within the user provided
    stream object and returns PSDs, PSD PDFs, and their corresponding statistics and masks (if desired).

    :param st: ObsPy Stream object
    :type st: obspy.core.stream.Stream
    :param expLoPeriod: Low end of the period band used for calculating the exponential dead channel detector
        (default = 4/trace sampling rate)
    :type expLoPeriod: int
    :param expHiPeriod: High end of the period band used for calculating the exponential dead channel detector
    :type expHiPeriod: int
    :param linLoPeriod: Low end of the period band used for calculating the linear dead channel detector
        (default = 4/trace sampling rate)
    :type linLoPeriod: int
    :param linHiPeriod: High end of the period band used for calculating the linear dead channel detector
    :type linHiPeriod: int
    :param evalresp: Response information obtained from numpy.ndarray of freq, amp, phase information
                     matching output of client.evalresp, optional for psdStatistics. If None, will pull from
                     `client.evalresp`
    :type evalresp: numpy.ndarray
    :param generateMasks: If True, boolean masks will be generated
    :type generateMasks: bool
    :param dcExpThreshold: Dead channel exponent threshold
    :type dcExpThreshold: int
    :param dcExpThresholdHour: Dead channel exponent hourly threshold. Only used if byHourOn is True
    :type dcExpThresholdHour: int
    :param pctBelowNoiseThreshold: Percent below NLNM threshold
    :type pctBelowNoiseThreshold: int
    :param pctAboveNoiseThreshold: Percent above NLNM threshold
    :type pctAboveNoiseThreshold: int
    :param pctBelowNoiseThresholdRESP: Percent below NLNM for bad response threshold
    :type pctBelowNoiseThresholdRESP: int
    :param pctAboveNoiseThresholdRESP: Percent above NLNM for bad response threshold.
    :type pctAboveNoiseThresholdRESP: int
    :param rmsThreshold: RMS threshold
    :type rmsThreshold: int
    :param dcLinThreshold: Dead channel linear threshold
    :type dcLinThreshold: int
    :param dcLinThresholdHour: Dead channel linear hourly threshold. Only used if byHouron set to True
    :type dcLinThresholdHour: int
    :param processes: Number of processes to use for calculation. 6 is the recommendation for personal desktop computer.
    :type processes: int
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger
    :param masksByTime: Boolean to determine whether masks are generated by time. If True, masks will be
                                    generated with a start/end time, if false, they will be generated as boolean array.
    :type masksByTime: bool
    :param byHourOn: Turn on hourly thresholding checks for dead_channel_exponential, dead_channel_linear,
                                 and dead_channel_gsn checks
    :type byHourOn: bool
    :param database_config: dictionary containing the necessary parameters to create
                            a pycheron Database object. 
                            These include "db_name", "session_name", "overwrite", "manual", "wfdb_conn"
    :type database_config: dict

    :return: List of dictionaries containing the following keys and types:

                   * start_time (`str`)
                   * end_time (`str`)
                   * snclq (`str`)
                   * percent_above_nhnm (`float`)
                   * percent_below_nlnm (`float`)
                   * dead_channel_exponent (`float`)
                   * dead_channel_linear (`float`)
                   * dead_channel_gsn (`bool`)
                   * uncorrected_psds (`list`)
                   * corrected_psds (`pandas.DataFrame`)
                   * pdfs (`pandas.DataFrame`)
                   * metric_name (`str`)
                   * dead_channel (`bool`)
                   * low_amp (`bool`)
                   * noise1(`bool`)
                   * noise2 (`bool`)
                   * highAmp (`bool`)
                   * badResp (`bool`)
                   * dc_mask (`numpy.ndarray`)
                   * dead_channel_gsn_mask (`numpy.ndarray`)
                   * low_amp_mask (`numpy.ndarray`)
                   * noise1_mask (`numpy.ndarray`)
                   * noise2_mask (`numpy.ndarray`)
                   * hi_amp_mask (`numpy.ndarray`)
                   * bad_resp_mask (`numpy.ndarray`)


    The following additional entries are included when byDayHourOn is `True`:
        * dead_channel_exponent_hourly (`list`)
        * dead_channel_linear_hourly (`list`)
        * dead_channel_gsn_hourly (`list`)
        * dead_chan_exp_hourly_masks (`numpy.ndarray`)
        * dead_chan_lin_hourly_masks (`numpy.ndarray`)
        * dead_chan_gsn_hourly_masks (`numpy.ndarray`)

    :rtype: list

    * Code originally ported from IRISMustangMetrics R Cran Package
    (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, March 20). CRAN-Package IRISMustangMetrics.
    The Comprehensive R Archive Network. Retrieved from
    https://cran.r-project.org/web/packages/IRISMustangMetrics/index.html) and augmented and adapted for use within
    Pycheron

    **Detailed Return Values:**

    .. list-table::
        :widths: 25 25 25 45
        :header-rows: 1

        * - Key
          - Value
          - Type
          - Details
        * - starttime
          - trace start time
          - `str`
          -
        * - endtime
          - trace end time
          - `str`
          -
        * - snclq
          - Station.Network.Channel.Location.Quality
          - `str`
          -
        * - pct_above_nhnm
          - percent above new high noise model
          - `float`
          - "Percentage of PSD values that are above the new High Noise model for their frequency.
            Only frequencies less than `sample_rate`/2 are considered to avoid instrument response effects as you
            approach the nyquist frequency. This value is calculated over the entire time period"
            (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)
        * - pct_below_nlnm
          - percent below new low noise model
          - `float`
          - "Percentage of PSD values that are below the new low noise model for their frequency.
             Only frequencies less than `sample_rate`/2 are considered to avoid instrument response
             effects as you approach the nyquist frequency.
             This value is calculated over the entire time period
            (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)
        * - dead_channel_exponent
          - dead channel metric exponential fit
          - `float`
          -  A "dead channel" metric is calculated from the mean of all the PSDs generated
                (typically 47 for a 24 hour period).
             Values of the PSD mean line over the band (expLoPeriod:expHiPeriod) are fit to an exponential.
             The `dead_channel_exp metric` is the standard deviation of the fit residuals.
             Lower numbers indicate a better fit and a higher likelihood that the mean PSD is exponential
                which is an indication of a "dead channel"
        * - dead_channel_exponent_hourly
          - dead channel metric exponential fit for each hour PSD
          - `list`
          -  A "dead channel" metric is calculated from each hourly PSD generated (typically 47 for a 24 hour
             period). Values of the PSD line over the band (expLoPeriod:expHiPeriod (defaults to 1 at high period band)
             are fit to an exponential. The `dead_channel_exp metric` is the standard deviation of the fit residuals
              (RMS error).
             Lower numbers indicate a better fit and a higher likelihood that the mean PSD is exponential
              -- an indication of a "dead channel"
        * - dead_channel_linear
          - dead channel metric linear fit
          - `float`
          - "A "dead channel" metric is calculated from the mean of all the PSDs generated (typically 47 for a 24 hour
            period). Values of the PSD mean line over the band (linLoPeriod:linHiPeriod) are fit to a line. The
            dead_channel_lin metric is the standard deviation of the fit residuals.
            Lower numbers indicate a better fit and a higher likelihood that the mean PSD is linear
                -- an indication of a "dead channel""
            (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)
        * - dead_channel_linear_hourly
          - dead channel metric linear fit for each hour PSD
          - `list`
          - "A "dead channel" metric is calculated from each" of hourly PSD generated "(typically 47 for a 24 hour
            period). Values of the PSD line over the band (linLoPeriod:linHiPeriod) are fit to a line. The
            dead_channel_lin metric is the standard deviation of the fit residuals (RMS error). Lower numbers indicate a
            better fit and a higher likelihood that the mean PSD is linear -- an indication of a "dead channel""
            (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)
        * - dead_channel_gsn
          - dead channel metric gsn fit
          - `bool`
          - Boolean metric using '0' or '1' to indicate if the average median deviation below the NLNM is sufficient to
            label channel as dead or not. If average > 5.0, mark dead_channel_gsn value = 1 (TRUE) (dead), else
            mark value = 0 (FALSE) (not dead).
            (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)
        * - dead_channel_gsn_hourly
          - dead channel metric gsn fit for each hour PSD
          - `list`
          - Boolean metric using '0' or '1' to indicate if the average median deviation below the NLNM is sufficient to
            label channel as dead or not for each hourly PSD. If average > 5.0, mark dead_channel_gsn value = 1 (TRUE)
            (dead), else mark value = 0 (FALSE) (not dead).
        * - uncorrected_psds
          - list of uncorrected PSDs for each trace in stream as output from the psdList function
          - `list`
          -
        * - corrected_psds
          - dataframe of starttime, endtime, frequency (Hz), power (dB) values of corrected PSDs
          - `pandas.DataFrame`
          -
        * - pdfs
          - dataframe of frequency (Hz), hits (count), power (dB), and power_mode (dB) values for PDFs.
            Power mode is used for the creation of colorgrids
          - `pandas.DataFrame`
          -
        * - metric_name
          - name of metric
          - `str`
          - psdMetric
          -
        * - dead_channel
          - Boolean value indicating whether dead channel threshold exceeded
            (corresponds to `dead_channel_exp` < `dcExpThreshold` and `avg_pctBelow` > `pctBelowNoiseThreshold`)
          - `bool`
          -
        * - low_amp
          - Boolean value indicating whether low amplitude threshold exceeded
            (corresponds to `dead_channel_lin` <= `dcLinThreshold` and `avg_pctBelow` > `pctBelowNoiseThreshold`)
          - `bool`
          -
        * - noise1
          - Boolean value indicating whether noise1 threshold exceeded
            (corresponds to `dead_channel_exp` < `dcExpThreshold` and `avg_pctAbove` > `pctAboveNoiseThreshold`)
          - `bool`
          -
        * - noise2
          - Boolean value indicating whether noise2 threshold exceeded
            (corresponds to `dead_channel_lin` >= `dcLinThreshold` and `avg_pctAbove` > `pctAboveNoiseThreshold`)
          - `bool`
          -
        * - highAmp
          - Boolean value indicating whether high amp threshold exceeded
            (corresponds to `rms` > `rmsThreshold` and `avg_pctAbove` > `pctAboveNoiseThreshold`)
          - `bool`
          -
        * - badResp
          - Boolean value indicating whether bad resp threshold exceeded
            (corresponds to `avg_pctAbove` > `pctAboveNoiseThresholdRESP` or `
             avg_pctBelow` > `pctBelowNoiseThresholdRESP`)
          - `bool`
          -
        * - dc_mask
          - dead channel mask
            (corresponds to `dead_channel_exp` < `dcExpThreshold` and `avg_pctBelow` > `pctBelowNoiseThreshold`)
          - `numpy.ndarray`
          -
        * - dead_channel_gsn_mask
          - dead channel gsn mask (corresponds to `dead_channel_gsn` == 1)
          - `numpy.ndarray`
          -
        * - low_amp_mask
          - low amplitude mask
            (corresponds to `dead_channel_lin` <= `dcLinThreshold` and `avg_pctBelow` > `pctBelowNoiseThreshold`)
          - `numpy.ndarray`
          -
        * - noise1_mask
          - noise mask 1
            (corresponds to `dead_channel_exp` < `dcExpThreshold` and `avg_pctAbove` > `pctAboveNoiseThreshold`)
          - `numpy.ndarray`
          -
        * - noise2_mask
          - noise mask 2
            (corresponds to `dead_channel_lin` >= `dcLinThreshold` and `avg_pctAbove` > `pctAboveNoiseThreshold`)
          - `numpy.ndarray`
          -
        * - hi_amp_mask
          - high amplitude mask (corresponds to `rms` > `rmsThreshold` and `avg_pctAbove` > `pctAboveNoiseThreshold`)
          - `numpy.ndarray`
          -
        * - bad_resp_mask
          - bad response mask (corresponds to `avg_pctAbove` > `pctAboveNoiseThresholdRESP` or
            `avg_pctBelow` > `pctBelowNoiseThresholdRESP`
          - `numpy.ndarray`
          -
        * - dead_chan_exp_hourly_masks
          - dead channel exponential hourly masks, masks generated when `dead_channel_exp_hou` <= `dcExpThresholdHour`
          - `numpy.ndarray`
          -
        * - dead_chan_lin_hourly_masks
          - dead channel linear hourly masks, masks generated when `dead_channel_lin_hour` <= `dcLinThresholdHour`
          - `numpy.ndarray`
          -
        * - dead_chan_gsn_hourly_masks
          - dead channel gsn hourly masks, masks generated when `dead_channel_gsn_hour` == 1
          - `numpy.ndarray`
          -

    **Algorithm Steps:**

    Algorithm steps below copied directly from IRISMustangMetrics code when ported from R to Python.

    #. Divide trace into Z-hour segments with 50% overlap

    #. For each Z-hour segment

       * Truncate segment to nearest power of 2 samples
       * Generate an averages/smoothed PSD as follows

         * Divide each truncated Z-hour segment into 13 segments with 75% overlap (Z/4 seconds length each)
         * For each of the 13 segments

           * Demean
           * Detrend
           * Apply 10% cosine taper
           * FFT
           * PSD
           * Normalize the power at each PSD frequency, multiplying it by (2*dt/Nseg), where Nseg is the # of
             samples in the segment

         * Average the 13 resulting PSDs to get averaged PSD
         * Multiply averaged PSD by 10% cosine taper scale factor (= 1.142857)
         * Frequency smooth the averaged PSD over 1-octave intervals at 1/8 octave increments (reducing # freq samples
           by a factor of 169; include 0.1 Hz as one of the geometric mean f values to sync f sampling)
         * Store the smoothed, averaged Z-hour PSD

    **Data Structure**

    .. code-block:: console

      [
        trace 1 in stream {
            start_time: starttime,
            end_time: endtime,
            snclq: snclq,
            percent_above_nhnm: avg_pctAbove,
            percent_below_nlnm: avg_pctBelow,
            dead_channel_exponent: dead_channel_exp,
            dead_channel_linear: dead_channel_lin,
            dead_channel_gsn: dead_channel_gsn,
            uncorrected_psds: [
                psd 1 [
                    [frequency array],
                    [psd array],
                    [snclq],
                    starttime,
                    endtime
                ],
                psd 2 [
                    [frequency array],
                    [psd array],
                    [snclq],
                    starttime,
                    endtime
                ],
                .
                .
                .
                psd 47 [
                    ...
                ]
            ],
            corrected_psds: DataFrame,
            pdfs: DataFrame,
            metric_name: "psdMetric",
            dc_mask: dcMask,
            low_amp_mask: lowAmpMask,
            noise1_mask: noise1Mask,
            noise2_mask: noise2Mask,
            hi_amp_mask: hiAmpMask,
            bad_resp_mask: badRESPMask
        },
        trace 2 in stream{
            ...
        }
      ]


    **Examples**

    .. code-block:: python

        #test data
        data = 'test/test_data/7a_cabn_bhe.884965.tar.mseed'
        #reading in stream
        st = obspy.read(data)
        #setting parameters
        expLoPeriod=None #Default
        expHiPeriod=100 #Default
        linLoPeriod=None #Default
        linHiPeriod=50 #Default
        evalresp=None
        #calculating metrics
        psd = psdMetric(st, expLoPeriod, expHiPeriod, linLoPeriod, linHiPeriod,evalresp)

        #print out example entries for the metrics,
        print 'Startime:', psd[0]['start_time']
        >>> Startime: 2013-11-01T00:00:00.000000Z
        print 'Endtime:', psd[0]['end_time']
        >>> Endtime: 2013-11-02T00:00:10.225000Z
        print 'SNCLQ:', psd[0]['snclq']
        >>> SNCLQ: 7A.CABN..BHE
        print 'Percent_above_nhnm:', psd[0]['percent_above_nhnm']
        >>> Percent_above_nhnm: 14.940378770166003
        print 'Percent_below_nlnm:', psd[0]['percent_below_nlnm']
        >>> Percent_above_nhnm: 14.940378770166003
        print 'Dead_Channel_Exponent:', psd[0]['dead_channel_exponent']
        >>> Dead_Channel_Exponent: 0.43778437537150205
        print 'Dead_Channel_Linear:', psd[0]['dead_channel_linear']
        >>> Dead_Channel_Linear: 5.880220527414604
        print 'Dead_Channel_GSN:', psd[0]['dead_channel_gsn']
        >>> Dead_Channel_GSN: 0
        print psd[0] # output not shown for brevity

    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Error check for empty stream objects
    if not st:
        logger.warn("psdMetric(): stopping PSD calculation because stream is empty")
        return

    # Initialize empty list to append information to
    d = []

    # Use psdList() function to apply the McNamara algorithm and compute some basic stats with the psdStatistics()
    # function
    # Note: All details about choosing window sizes, etc. are done in the psdList function. See psdList and McNamara
    # PSD method for more details on algorithm
    psds = psdList(st, processes=processes, logger=logger)
    # Let users know psdList calculated
    logger.log("psdMetric(): psdList Completed")

    # if psds exists, check if evalresp provided. If it hasn't, psdStatistics will go out to IRIS and fetch the info it
    # needs, else, use the user provided response
    if psds:
        if evalresp is None:
            stats = psdStatistics(psds, evalresp=None, logger=logger)
        else:
            stats = psdStatistics(psds, evalresp=evalresp, logger=logger)
    # If nothing returned from psdList, error out
    else:
        logger.error("psdMetric(): No PSDS Computed")
        return

    # Check if stats is empty, if so exit
    if not stats:
        logger.error("psdMetric(): Stats empty, exiting")
        return

    # get rms for high amp rms masks
    basicStats = basicStatsMetric(st, rmsThreshold=rmsThreshold)
    rms = basicStats[0]["rms"]

    # Loop through psd stats and grab out useful metadata information, such as starttime, endtime, sampling rate.
    # Essentially, if there is more than one snclq in the stream object this will loop through the various corresponding
    # stats/psds and utilize the appropriate ones. Also calculates the expLo/linLoPeriods based on the appropriate
    # streams
    for i in range(len(stats)):
        # initialize masks
        dcMask = None
        deadChanGSNMask = None
        lowAmpMask = None
        noise1Mask = None
        noise2Mask = None
        hiAmpMask = None
        badRESPMask = None
        dcExpByHourMask = None
        dcLinByHourMask = None
        dcGSNByHourMask = None

        # Get trace info, starttime, endtime, snclq, sampling rate
        tr = st[i]
        starttime = tr.stats.starttime
        endtime = tr.stats.endtime
        snclq_tr = st[i].get_id()
        trStats = stats[i]
        un_psds = psds[i]
        samplingRate = tr.stats.sampling_rate
        # Calculate expLoPeriod and linLoPeriod based on sampling rate
        if expLoPeriod is None:
            expLoPeriod = 4 / samplingRate
        if linLoPeriod is None:
            linLoPeriod = 4 / samplingRate

        # "pctAbove, pctBelow are functions of frequency returned by psdStats
        # Only frequencies less than nyquist/1.5 are used when calculating avg_pctAbove and avg_pctBelow
        # due to instrument response effects on PSD power as you approach the nyquist (= sample rate/2)
        # This is reasonable for most modern use cases, excluding older seismometers that may
        # have relied upon anti-alias filters. We do not seek to accommodate these at this point in time."
        # (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)

        nyquist = samplingRate / 2
        avg_pctAbove = np.nanmean(trStats["percent_above_nhnm"][trStats["noise_matrix_frequency"] < nyquist / 1.5])
        avg_pctBelow = np.nanmean(trStats["percent_below_nlnm"][trStats["noise_matrix_frequency"] < nyquist / 1.5])

        # Dead Channel Exponential fit metric
        # "The dead_channel_exponential metric is calculated by fitting the PSD mean line as seen in a PDF plot to an
        # exponential and calculating the standard deviation of the residuals.
        # The mean of a healthy set of PSDs will have a very non-exponential shape and large residuals while a
        # "dead channel" (aka digitizer noise) will have a PSD mean that appears as an exponential decay as a function
        # of log10(period).
        # The dead_channel_exp metric looks for an overall exponential decay regardless of the frequency range of the
        # detector. As the frequency bands will be sensor specific because of the calculation in the psdList() function,
        # we can by with just lopping an integer number of bands rather than specifying a band pass region we want.
        # IRIS says this algorithm is purely heuristic and resulted from a visual assessment of PDF plots of channels
        # known to be "dead".
        # The dead_channel_exp and dead_channel_lin metrics are not valid for LH? or VH? channels
        #
        # Exp fit steps:
        #   1) Determine index range (inverted because of freq->period conversion)
        #   2) Convert band from PSD mean to positive, non-zero values
        #   3) Fit log10(PSD mean) vs. log10(period) to a line
        #   4) Calculate the standard deviation of residuals of the fit(How close to exponential is the PSD mean line?)"
        #   (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)

        # Convert frequency to period
        period = 1 / trStats["noise_matrix_frequency"]

        # Determine index range of period values that are greater than or equal to expHiPeriod and less than or equal
        # to expLoPeriod.
        try:
            first = np.max([i for i, v in np.ndenumerate(period) if v >= expHiPeriod]) + 1
            if byHourOn:
                # This cuts off the latter half of the PSD, which seems to work better, especially if there is concavity
                # expHiPeriod = 1 //TODO this is causing errors where the first index is greater than the last
                first = np.max([i for i, v in np.ndenumerate(period) if v >= expHiPeriod]) + 1
        # Warn users if Max value could not be computed
        except ValueError:
            logger.warn("psdMetric(): Value Error - Max value could not be computed.")
            continue
        try:
            last = np.min([i for i, v in np.ndenumerate(period) if v <= expLoPeriod])
        # Warn users if min value could not be computed
        except ValueError:
            logger.warn("psdMetric(): Value Error - Min value could not be computed.")
            continue

        # If doing deadChanExponential by hours too, to try to get more granularity on issues
        if byHourOn:
            dead_channel_exp_hour = []
            for j in range(len(stats[0]["noise_matrix_noise"])):
                # Adjust PSD to positive values for log calculation, if values still less than one, take min again.
                # Subtracting from min of mean instead of min of itself is because for the PSD ones tested subtracting
                # from itself seemed to actually change the shape of the PSD, which is contrary to what we are trying to
                # do
                try:
                    positivePSD = stats[0]["noise_matrix_noise"][j][first:last] - min(trStats["mean"][first:last])
                    if np.where(positivePSD < 1):
                        positivePSD = positivePSD - min(positivePSD) + 0.1

                    # Exponential fit; fit log10(PSD) vs. log10(period) to a line; used polyfit here instead of
                    # np.linalg.lstsq as R codes are simply fitting to line,
                    # and polyfit will take the input directly in 1D.
                    # np.linalg.lstsq creates more output than is necessary (condenses code from 8 lines to 1 line)
                    a = np.polyfit(np.log10(period[first:last]), np.log10(positivePSD), 1)

                    # Calculate standard deviation of the residuals of the fits to obtain dead_channel_exp metric
                    dead_channel_exp_hour.append(
                        np.std((np.log10(positivePSD) - (np.log10(period[first:last]) * a[0] + a[1])))
                    )

                    # Exp fit is not valid for data 1sps and below, in order to ensure this isn't used in
                    # thresholding/throwing out data for such channels we do the lazy way and simply overwrite it as
                    # 'Not Valid', as data >= 1 sps can still be used with GSN metric,
                    # thus we don't want to limit the 1sps
                    # and less data from continuing on with the subsequent calculations
                    if samplingRate <= 1:
                        dead_channel_exp_hour.append("Metric not valid for sample rate")
                except (ValueError, TypeError):
                    continue

        # Convert band from PSD mean to positive, non-zero values
        positiveMean = trStats["mean"][first:last] - min(trStats["mean"][first:last]) + 0.1

        # Exponential fit; fit log10(PSDmean) vs. log10(period) to a line; used polyfit here instead of np.linalg.lstsq
        # as R codes are simply fitting to line, and polyfit will take the input directly in 1D. np.linalg.lstsq creates
        # more output than is necessary (condenses code from 8 lines to 1 line)
        a = np.polyfit(np.log10(period[first:last]), np.log10(positiveMean), 1)

        # Calculate standard deviation of the residuals of the fits (RMS error) to obtain dead_channel_exp metric
        dead_channel_exp = np.std((np.log10(positiveMean) - (np.log10(period[first:last]) * a[0] + a[1])))

        # Exp fit is not valid for data 1sps and below, in order to ensure this isn't used in thresholding/throwing out
        # data for such channels we do the lazy way and simply overwrite it as 'Not Valid', as data >= 1 sps can still
        # be used with GSN metric, thus we don't want to limit the 1sps and less data from continuing on with the
        # subsequent calculations
        if samplingRate <= 1:
            dead_channel_exp = "Metric not valid for sample rate"

        # "Linear fit: Fits the PSD mean line as a linear function of log10(period) over a specific
        # band. Some types of sensor malfunction result in a flat spectra that is different
        # from a curve produced by digitizer noise. Be aware that some normal noise patterns will also score low on this
        # metric.
        #
        # Linear fit steps:
        #   1) Determine index range (inverted because of freq->period conversion)
        #   2) Fit PSD mean vs. log10(period) to a line
        #   3) Calculate the standard deviation of residuals of the fit
        # The dead_channel_exp and dead_channel_lin metrics are not valid for LH? or VH? channels"
        # (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)

        # Determine index range of period values that are greater than or equal to linHiPeriod and less than or equal
        # to linLoPeriod.
        try:
            firstLin = np.max([i for i, v in np.ndenumerate(period) if v >= linHiPeriod]) + 1
        # Warn users if max value could not be computed
        except ValueError:
            logger.warn("psdMetric(): Value Error - Max value could not be computed.")
            continue
        try:
            lastLin = np.min([i for i, v in np.ndenumerate(period) if v <= linLoPeriod])
        # Warn users if min value could not be computed
        except ValueError:
            logger.warn("psdMetric(): Value Error - Min value could not be computed.")
            continue

        # If doing deadChanLinear by hours too, to try to get more granularity on issues
        if byHourOn:
            dead_channel_lin_hour = []
            for j in range(len(stats[0]["noise_matrix_noise"])):
                try:
                    # Grab out values using firstlin and lastlin for each PSD
                    psdLin = stats[0]["noise_matrix_noise"][j][firstLin:lastLin]

                    # Linear fit; fit psd vs. log10(period) to a line. used polyfit here instead of np.linalg.lstsq
                    # as R codes are simply fitting to line, and polyfit will take the input directly in 1D.
                    # np.linalg.lstsq creates more output than is necessary (condenses code from 3 lines to 1 lines)
                    a = np.polyfit(np.log10(period[firstLin:lastLin]), psdLin, 1)

                    # Calculate standard deviation of the residuals of the fits (RMS error) to obtain dead_channel_lin
                    # metric
                    dead_channel_lin_hour.append(np.std(psdLin - (np.log10(period[firstLin:lastLin]) * a[0] + a[1])))

                    # Lin fit is not valid for data 1sps and below, in order to ensure this isn't used in
                    # thresholding/throwing out data for such channels we do the lazy way and simply overwrite it as
                    # 'Not Valid', as data >= 1 sps can still be used with GSN metric,
                    # thus we don't want to limit the 1sps
                    # and less data from continuing on with the subsequent calculations
                    if samplingRate <= 1:
                        dead_channel_lin_hour.append("Metric not valid for sample rate")
                except (ValueError, TypeError):
                    continue

        # Grab out appropriate psdMean values using firstlin and lastLin
        psdMean = trStats["mean"][firstLin:lastLin]

        # Linear fit; fit psd mean vs. log10(period) to a line. used polyfit here instead of np.linalg.lstsq
        # as R codes are simply fitting to line, and polyfit will take the input directly in 1D. np.linalg.lstsq creates
        # more output than is necessary (condenses code from 3 lines to 1 lines)
        a = np.polyfit(np.log10(period[firstLin:lastLin]), psdMean, 1)

        # Calculate standard deviation of the residuals of the fits (RMS error) to obtain dead_channel_lin metric
        dead_channel_lin = np.std(psdMean - (np.log10(period[firstLin:lastLin]) * a[0] + a[1]))

        # Lin fit is not valid for data 1sps and below, in order to ensure this isn't used in thresholding and throwing
        # out data for such channels, we do the lazy way and simply overwrite it as 'Not Valid', as data >= 1sps can
        # still be used with GSN metric, thus we don't want to limit the 1sps and less data from continuing on with
        # calculations
        if samplingRate <= 1:
            dead_channel_lin = "Metric not valid for sample rate"

        # Note that the dead_channel_exp and dead_channel_lin metrics are not valid for LH? or VH? channels,
        # thus another metric was added to deal with looking at a narrow band of the PDF to determine if
        # the deviation below the NLNM is sufficient to label it as dead or not dead. This is a boolean metric using
        # numeric indicators, called dead_channel_gsn
        #
        # "GSN dead channel metric steps:
        #  1) Only accept 1 SPS and greater sample rate
        #  2) Get the PDF matrix and slice the 4 to 8 second period
        #  3) Calculate the median value at each period step, writing result to a vector
        #  4) Collect differences in the median vector from the equivalent NLNM vector [ NLNM(x) - Y(x) ]
        #  5) Average the result of the differences
        #  6) If average > 5.0, mark dead_channel_gsn value = 1 (TRUE), else mark value = 0 (FALSE)"
        # (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)

        dead_channel_gsn = 0

        # Use condition that sample rate must be >= 1 sps
        if samplingRate >= 1:
            # Grab out values of the PDF noise matrix >= 4 seconds and <= 8  seconds period
            firstGSN = np.max([i for i, v in np.ndenumerate(period) if v >= 4]) + 1
            firstGSN += 1
            lastGSN = np.min([i for i, v in np.ndenumerate(period) if v <= 8]) - 1
            psdSlice = trStats["noise_matrix_noise"][:, firstGSN:lastGSN:-1]
            unH_floor = np.floor(np.min(psdSlice))
            pdfSlice = trStats["pdf_matrix"][:, firstGSN:lastGSN:-1]

            # Get column-wise median of 'counted' binned dB values
            psdMedian = np.median(unHistogram(pdfSlice, unH_floor, 1), axis=0)

            # Ensure we have an average low dB value and if so, treat as a dead channel condition
            # Grab out NLNM slice to match, reverse order to match how pulled in R code base
            nlnmSlice = trStats["nlnm"][lastGSN:firstGSN][::-1]

            # Get difference between NLNM and values in psdMedian, then get average of differences (i.e., deviation)
            diffToNM = nlnmSlice[~np.isnan(psdMedian)] - psdMedian
            avgDiff = np.mean(diffToNM)

            # Compare average to 5 db, save 1 if condition is true, and 0 if condition is false
            if avgDiff > 5:
                dead_channel_gsn = 1
            else:
                dead_channel_gsn = 0

            # If byHourOn, loop through individual PSD slices
            if byHourOn:
                # Create empty lists to append to
                dead_channel_gsn_hour = []
                avgDiff = []
                # Loop through individual PSD slices for the period band of interest pulled above
                for j in range(len(psdSlice)):
                    # Difference PSD from nlnmSlice for same periods, then take average
                    diffToNM = nlnmSlice - psdSlice[j]
                    avgDiff.append(np.mean(diffToNM))

                # Loop through average diff and put
                for j in range(len(avgDiff)):
                    if avgDiff[j] > 5:
                        dead_channel_gsn_hour.append(1)
                    else:
                        dead_channel_gsn_hour.append(0)

            # Create a corrected PSD data frame similar to that returned by
            # http://service.iris.edu/mustang/noise-psd/1/
            # Grab out all start times and endtimes from psds (output from psdList)
            startlist = []
            endlist = []
            for j in range(len(psds[0])):
                startlist.append(psds[0][j][3])
                endlist.append(psds[0][j][4])
            # Create repeating elements of frequency based on number of rows in trStats['noise_matrix_noise']
            # (output from PSDStatistics)
            freqV = np.repeat(trStats["noise_matrix_frequency"], len(trStats["noise_matrix_noise"]))

            # Concatenate power values from trStats['noise_matrix_noise'] (output from PSDStatistics) into single array
            powerV = np.concatenate(trStats["noise_matrix_noise"])

            # Create list of repeated values of start list based on the number of columns in noise matrix and length of
            # startlist
            startV = np.repeat(
                startlist,
                np.repeat(len(trStats["noise_matrix_noise"][0]), len(startlist)),
            )
            endV = np.repeat(endlist, np.repeat(len(trStats["noise_matrix_noise"][0]), len(endlist)))

            # Finally, now create the dataframe of corrected PSDs
            cPsdDF = pd.DataFrame(
                OrderedDict(
                    (
                        ("starttime", pd.Series(startV)),
                        (("endtime", pd.Series(endV))),
                        (("frequency", pd.Series(freqV))),
                        (("power", pd.Series(powerV))),
                    )
                )
            )

            # Create similar dataframe for PDFs
            # Concatenate power values from trStats['PDFMatrix'] (output from PSDStatistics) into single array
            hitsV = np.concatenate(trStats["pdf_matrix"])

            # Create repeating elements of PDFBins based on number of rows in trStats['PDFMatrix']
            # (output from PSDStatistics)
            binV = np.repeat(trStats["pdf_bins"], len(trStats["pdf_matrix"][0]))

            # Create repeating elements of frequency based on number of rows in trStats['PDFMatrix']
            # (output from PSDStatistics) and length of trStats['noise_matrix_frequency']
            freq2V = np.repeat(
                trStats["noise_matrix_frequency"],
                np.repeat(len(trStats["pdf_matrix"]), len(trStats["noise_matrix_frequency"])),
            )

            # Create repeating elements of power mode based on number of rows in trStats['PDFMatrix']
            # (output from PSDStatistics) and length of trStats['Mode']
            modeV = np.repeat(
                trStats["mode"],
                np.repeat(len(trStats["pdf_matrix"]), len(trStats["noise_matrix_frequency"])),
            )

            # Create dataframe now
            pdfDF = pd.DataFrame(
                data={
                    "frequency": freq2V,
                    "power": binV,
                    "hits": hitsV,
                    "power_mode": modeV,
                }
            )

            # Filter out values in the dataframe to only include values where hits > 0
            pdfDf = pdfDF[pdfDF["hits"] > 0]

            pctBelowNLNM = np.repeat(trStats["percent_below_nlnm"], len(trStats["noise_matrix_noise"]))
            pctAboveNLNM = np.repeat(trStats["percent_above_nhnm"], len(trStats["noise_matrix_noise"]))

            # Boolean threshold checks
            # dead channel check
            if dead_channel_exp < dcExpThreshold and avg_pctBelow > pctBelowNoiseThreshold:
                dead_channel = 1
            else:
                dead_channel = 0
            # low amp check
            if dead_channel_lin <= dcLinThreshold and avg_pctBelow > pctBelowNoiseThreshold:
                low_amp = 1
            else:
                low_amp = 0
            # noise 1
            if dead_channel_exp < dcExpThreshold and avg_pctAbove > pctAboveNoiseThreshold:
                noise1 = 1
            else:
                noise1 = 0
            # noise 2
            if dead_channel_lin >= dcLinThreshold and avg_pctAbove > pctAboveNoiseThreshold:
                noise2 = 1
            else:
                noise2 = 0
            # high amp
            if rms > rmsThreshold and avg_pctAbove > pctAboveNoiseThreshold:
                highAmp = 1
            else:
                highAmp = 0
            # bad resp
            if avg_pctAbove > pctAboveNoiseThresholdRESP or avg_pctBelow > pctBelowNoiseThresholdRESP:
                badResp = 1
            else:
                badResp = 0

            if generateMasks:

                if not masksByTime:
                    # dead channel mask
                    # Check threshold, if exceeds apply 1 to the entire thing
                    if dead_channel_exp < dcExpThreshold and avg_pctBelow > pctBelowNoiseThreshold:
                        # entire mask is 1
                        dcMask = np.ones(len(cPsdDF))

                    # dead channel gsn mask
                    # If dead channel gsn triggered, apply 1 to entire thing
                    if dead_channel_gsn == 1:
                        deadChanGSNMask = np.ones(len(cPsdDF))

                    # low Amp
                    # Check threshold, if met, then apply 1 to indexes where low amp is occurring
                    if dead_channel_lin <= dcLinThreshold and avg_pctBelow > pctBelowNoiseThreshold:
                        lowAmpMask = np.zeros(len(cPsdDF))
                        index_pct = np.where(pctBelowNLNM > pctBelowNoiseThreshold)
                        lowAmpMask[index_pct] = 1

                    # noise1
                    # Check threshold, if met, then apply 1 to indexes where noise is occurring
                    if dead_channel_exp < dcExpThreshold and avg_pctAbove > pctAboveNoiseThreshold:
                        noise1Mask = np.zeros(len(cPsdDF))
                        index_pct = np.where(pctAboveNLNM > pctAboveNoiseThreshold)
                        noise1Mask[index_pct] = 1

                    # noise 2
                    # Check threshold, if met, then apply 1 to indexes where noise2 is occurring
                    if dead_channel_lin >= dcLinThreshold and avg_pctAbove > pctAboveNoiseThreshold:
                        noise2Mask = np.zeros(len(cPsdDF))
                        index_pct = np.where(pctBelowNLNM <= pctBelowNoiseThreshold)
                        noise2Mask[index_pct] = 1

                    # hiAmp
                    # Check threshold, if met, then apply 1 to indexes where high amp is occurring
                    if rms > rmsThreshold and avg_pctAbove > pctAboveNoiseThreshold:
                        hiAmpMask = np.zeros(len(cPsdDF))
                        index_pct = np.where(pctBelowNLNM > pctBelowNoiseThreshold)
                        hiAmpMask[index_pct] = 1

                    # bad resp
                    # Check threshold, if met, then apply 1 to indexes where bad resp is occurring
                    if avg_pctAbove > pctAboveNoiseThresholdRESP or avg_pctBelow > pctBelowNoiseThresholdRESP:
                        badRESPMask = np.zeros(len(cPsdDF))
                        index_pct = np.where(pctBelowNLNM > pctBelowNoiseThresholdRESP)
                        index_pct_1 = np.where(pctAboveNLNM > pctAboveNoiseThresholdRESP)
                        badRESPMask[index_pct] = 1
                        badRESPMask[index_pct_1] = 1

                    if byHourOn:  # This doesn't seem right
                        # Dead Channel Exponential masks by hour
                        dcExpByHourMask = np.zeros(len(cPsdDF))
                        deadChanExpHour = np.repeat(
                            dead_channel_exp_hour,
                            len(trStats["noise_matrix_frequency"]),
                        )
                        for j in range(len(deadChanExpHour)):
                            if deadChanExpHour[j] <= dcExpThresholdHour:
                                dcExpByHourMask[j] = 1
                            else:
                                dcExpByHourMask[j] = 0

                        # Dead Channel Linear masks by hour
                        dcLinByHourMask = np.zeros(len(cPsdDF))
                        deadChanLinHour = np.repeat(
                            dead_channel_lin_hour,
                            len(trStats["noise_matrix_frequency"]),
                        )
                        for j in range(len(deadChanLinHour)):
                            if deadChanLinHour[j] >= dcLinThresholdHour:
                                dcLinByHourMask[j] = 1
                            else:
                                dcLinByHourMask[j] = 0

                        # Dead Channel GSN masks by hour
                        dcGSNByHourMask = np.repeat(
                            dead_channel_gsn_hour,
                            len(trStats["noise_matrix_frequency"]),
                        )

                else:
                    # dead channel mask
                    if dead_channel_exp < dcExpThreshold and avg_pctBelow > pctBelowNoiseThreshold:
                        # entire mask is 1, if exceeded
                        dcMask = [
                            {
                                "starttime": starttime.isoformat(),
                                "endtime": endtime.isoformat(),
                            }
                        ]

                    # dead channel gsn mask
                    if dead_channel_gsn == 1:
                        # entire mask is 1 if boolean value = 1
                        deadChanGSNMask = [
                            {
                                "starttime": starttime.isoformat(),
                                "endtime": endtime.isoformat(),
                            }
                        ]

                    # low Amp
                    # Check threshold, if met, then create masks for times where issues occur
                    if dead_channel_lin <= dcLinThreshold and avg_pctBelow > pctBelowNoiseThreshold:
                        index_pct = np.where(pctBelowNLNM > pctBelowNoiseThreshold)[0]
                        lowAmpMask = samples2time(index_pct, samplingRate, starttime)

                    # noise1
                    # Check threshold, if met, then create masks for times where issues occur
                    if dead_channel_exp < dcExpThreshold and avg_pctAbove > pctAboveNoiseThreshold:
                        index_pct = np.where(pctAboveNLNM > pctAboveNoiseThreshold)[0]
                        noise1Mask = samples2time(index_pct, samplingRate, starttime)

                    # noise 2
                    # Check threshold, if met, then create masks for times where issues occur
                    if dead_channel_lin >= dcLinThreshold and avg_pctAbove > pctAboveNoiseThreshold:
                        index_pct = np.where(pctBelowNLNM <= pctBelowNoiseThreshold)[0]
                        noise2Mask = samples2time(index_pct, samplingRate, starttime)

                    # hiAmp
                    # Check threshold, if met, then create masks for times where issues occur
                    if rms > rmsThreshold and avg_pctAbove > pctAboveNoiseThreshold:
                        index_pct = np.where(pctBelowNLNM > pctBelowNoiseThreshold)[0]
                        hiAmpMask = samples2time(index_pct, samplingRate, starttime)

                    # bad resp
                    # Check threshold, if met, then create masks for times where issues occur
                    if avg_pctAbove > pctAboveNoiseThresholdRESP or avg_pctBelow > pctBelowNoiseThresholdRESP:
                        index_pct = np.where(pctBelowNLNM > pctBelowNoiseThresholdRESP)[0]
                        index_pct_1 = np.where(pctAboveNLNM > pctAboveNoiseThresholdRESP)[0]
                        badRESPMask = samples2time(index_pct, samplingRate, starttime)
                        badRESPMask = badRESPMask + samples2time(index_pct_1, samplingRate, starttime)

                    if byHourOn:
                        # Dead Channel Exponential masks by hour
                        deadChanExpHour = np.repeat(
                            dead_channel_exp_hour,
                            len(trStats["noise_matrix_frequency"]),
                        )
                        index_pct = np.where(deadChanExpHour <= dcExpThresholdHour)[0]
                        dcExpByHourMask = samples2time(index_pct, samplingRate, starttime)
                        if not dcExpByHourMask:
                            dcExpByHourMask = None

                        # Dead Channel Linear masks by hour
                        deadChanLinHour = np.repeat(
                            dead_channel_lin_hour,
                            len(trStats["noise_matrix_frequency"]),
                        )
                        index_pct = np.where(deadChanLinHour <= dcLinThresholdHour)[0]
                        dcLinByHourMask = samples2time(index_pct, samplingRate, starttime)
                        if not dcLinByHourMask:
                            dcLinByHourMask = None

                        # Dead Channel GSN
                        deadChanGSNHour = np.repeat(
                            dead_channel_gsn_hour,
                            len(trStats["noise_matrix_frequency"]),
                        )
                        index_pct = np.where(deadChanGSNHour == 1)[0]
                        dcGSNByHourMask = samples2time(index_pct, samplingRate, starttime)
                        if not dcGSNByHourMask:
                            dcGSNByHourMask = None
        # Define metrics dictionary if byHour on
        if byHourOn:
            metrics = {
                "start_time": starttime.isoformat(),
                "end_time": endtime.isoformat(),
                "snclq": snclq_tr,
                "percent_above_nhnm": avg_pctAbove,
                "percent_below_nlnm": avg_pctBelow,
                "dead_channel_exponent": dead_channel_exp,
                "dead_channel_exponent_hourly": dead_channel_exp_hour,
                "dead_channel_linear": dead_channel_lin,
                "dead_channel_linear_hourly": dead_channel_lin_hour,
                "dead_channel_gsn": dead_channel_gsn,
                "dead_channel_gsn_hourly": dead_channel_gsn_hour,
                "uncorrected_psds": un_psds,
                "corrected_psds": cPsdDF,
                "pdfs": pdfDf,
                "metric_name": "psdMetric",
                "dead_channel": dead_channel,
                "low_amp": low_amp,
                "noise1": noise1,
                "noise2": noise2,
                "highAmp": highAmp,
                "badResp": badResp,
                "dc_mask": dcMask,
                "dead_channel_gsn_mask": deadChanGSNMask,
                "low_amp_mask": lowAmpMask,
                "noise1_mask": noise1Mask,
                "noise2_mask": noise2Mask,
                "hi_amp_mask": hiAmpMask,
                "bad_resp_mask": badRESPMask,
                "dead_chan_exp_hourly_masks": dcExpByHourMask,
                "dead_chan_lin_hourly_masks": dcLinByHourMask,
                "dead_chan_gsn_hourly_masks": dcGSNByHourMask,
            }
        # Otherwise define metrics
        else:

            metrics = {
                "start_time": starttime.isoformat(),
                "end_time": endtime.isoformat(),
                "snclq": snclq_tr,
                "percent_above_nhnm": avg_pctAbove,
                "percent_below_nlnm": avg_pctBelow,
                "dead_channel_exponent": dead_channel_exp,
                "dead_channel_linear": dead_channel_lin,
                "dead_channel_gsn": dead_channel_gsn,
                "uncorrected_psds": un_psds,
                "corrected_psds": cPsdDF,
                "pdfs": pdfDf,
                "metric_name": "psdMetric",
                "dead_channel": dead_channel,
                "low_amp": low_amp,
                "noise1": noise1,
                "noise2": noise2,
                "highAmp": highAmp,
                "badResp": badResp,
                "dc_mask": dcMask,
                "dead_channel_gsn_mask": deadChanGSNMask,
                "low_amp_mask": lowAmpMask,
                "noise1_mask": noise1Mask,
                "noise2_mask": noise2Mask,
                "hi_amp_mask": hiAmpMask,
                "bad_resp_mask": badRESPMask,
            }
        # Append metrics to overall list object
        d.append(metrics)

    # If database defined, insert metric information
    if database_config is not None:
        database = Database(**database_config)
        database.insert_metric(d)

    return d
