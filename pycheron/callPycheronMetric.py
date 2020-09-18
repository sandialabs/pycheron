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

import warnings

warnings.filterwarnings("ignore")

import fnmatch
import os
import sys
import matplotlib.pyplot as plt
import time

import numpy as np
import obspy
from obspy.core.utcdatetime import UTCDateTime
from obspy.imaging.cm import pqlx

from pycheron.plotting.dailyPdfPlot import dailyPdfplots
from pycheron.plotting.psdPlot import psdPlot
from pycheron.plotting.rankingPlot import rankPlot
from pycheron.psd.noise.stationNoiseModel import stationNoiseModel
from pycheron.util.getChannel import *
from pycheron.util.getSta import *
from pycheron.util.getNet import *
from pycheron.util.format import parse_snclq
from pycheron.psd.noise.networkNoiseModel import networkNoiseModel
from pycheron.dataAcq.css import *
from pycheron.metrics.basicStatsMetric import basicStatsMetric
from pycheron.metrics.correlationMetric import correlationMetric
from pycheron.metrics.crossCorrMetric import crossCorrMetric
from pycheron.metrics.DCOffSetTimesMetric import DCOffSetTimesMetric
from pycheron.metrics.gapMetric import gapMetric
from pycheron.metrics.repeatedAmplitude import repeatedAmplitudeMetric
from pycheron.metrics.deadChannelMetric import deadChannelMetric
from pycheron.metrics.deadChanADFMetric import deadChanADFMetric
from pycheron.metrics.deadChanMeanMetric import deadChanMean
from pycheron.metrics.snrMetric import snrMetric
from pycheron.metrics.psdMetric import psdMetric
from pycheron.metrics.sohMetric import sohMetric
from pycheron.metrics.spikesMetric import spikesMetric
from pycheron.metrics.staltaMetric import staltaMetric
from pycheron.metrics.calibration import *
from pycheron.metrics.transferFunctionMetric import transferFunctionMetric
from pycheron.util.logger import Logger
from pycheron.db.sqllite_db import Database

from multiprocessing import Process


# -----------------------------------------------------------------------------
# Pycheron Code
# -----------------------------------------------------------------------------


def callPycheron(
    output_dir,
    data,
    datatype="dir",
    calcAll=True,
    calcPsds=False,
    calcBasic=False,
    calcCorr=False,
    calcCrossCorr=False,
    calcGap=False,
    calcAmp=False,
    calcSNR=False,
    calcSOH=False,
    calcStalta=False,
    calcDcOffset=False,
    calcSpikes=False,
    calcAllDeadChan=False,
    calcTransfer=False,
    calcCal=False,
    network=None,
    station=None,
    byDay=True,
    startdate=None,
    enddate=None,
    jul_start=None,
    jul_end=None,
    generateMasks=False,
    masksByTime=True,
    rmsThreshold=50000,
    maxThreshold=None,
    minThreshold=None,
    medianThreshold=None,
    meanThreshold=None,
    varianceThreshold=None,
    stdThreshold=None,
    maxLagSecs=10,
    filt="lowpass",
    freqmin=1.0,
    freqmax=100.0,
    corners=4,
    zerophase=False,
    maxorder=None,
    ba=False,
    freq_passband=False,
    windowSecs=1800,
    incrementSecs=None,
    threshold=0.9,
    separateMasks=True,
    completeDay=True,
    expLoPeriod=None,
    expHiPeriod=100,
    linLoPeriod=None,
    linHiPeriod=50,
    evalresp=None,
    dcExpThreshold=0.3,
    pctBelowNoiseThreshold=20,
    pctAboveNoiseThreshold=20,
    dcLinThreshold=2,
    num_gaps=10,
    pctBelowNoiseThresholdRESP=90,
    pctAboveNoiseThresholdRESP=90,
    minRep=10,
    algorithmSNR="splitWindow",
    windowSecsSNR=60,
    snrThreshold=2,
    data_quality=False,
    activity=False,
    io_clock=False,
    windowSize=41,
    thresholdSpikes=10,
    selectivity=None,
    fixedThreshold=True,
    staSecs=3,
    ltaSecs=30,
    increment=1,
    algorithmSTA="classic_LR",
    plots=True,
    pdfModel="gsn",
    per_arr=None,
    showNoiseModel=True,
    showMaxMin=False,
    showMode=True,
    showMean=False,
    showMedian=False,
    showEnvelope=False,
    envelopeType="10_90",
    showSingle=False,
    singleType=None,
    ylo=-200,
    yhi=-50,
    min_stations=5,
    rank_by=None,
    processesPSD=7,
    processesSpikes=7,
    log="pycheron.log",
    fortran=True,
    timespan="all",
    dcADF_win_size=500,
    dcADF_pval_thresh=0.01,
    dcADF_threshold=1.5,
    dcADF_use_thresh="pvalue",
    dcMean_win_size=None,
    dcMean_thresh=0.05,
    cal_metric_store=None,
    dcExpThresholdHour=0.25,
    dcLinThresholdHour=2,
    byHourOn=False,
    database="pycheron.db",
    session="None",
    overwrite=True,
    to_csv=False,
    stationStartAt=None,
):
    """
    Pycheron wrapper script that runs the entire program. See the :doc:`../tutorials` for an in
    depth tutorial of usage.

    **General**

    :param output_dir: Path to location where top level <network_dir> will be saved
    :type output_dir: string
    :param data: Data to be processed. Currently supports:

                * MSEED: (str) Path to single .mseed file
                * Wfdisc: (str) Path to .wfdisc file
                * Stream: Obspy Stream Object
                * Directory: (str) Path to directory of multiple .mseed files. With filename format:
                             ``<station>_<channel>_<julian day>.mseed``. Use :doc:`../pycheron.util.format` to format
                             data correctly.

    :param datatype: The data type of parameter *data*.

                    * "mseed": Use if *data* is single .mseed file
                    * "stream": Use if *data* is Obspy Stream Object
                    * "wfdisc": Use if *data* is path to .wfdisc file
                    * "dir": Use if *data* is path to directory.

    :type datatype: str
    :param calcAll: Calculate all metrics. If set to False, individual metrics will be ran based on
                           parameters ``calc<Metric>``.
    :type calcAll: bool
    :param network: Network of stations. **Note:** Parameter required when using wfdic format.
    :type network: str
    :param station: Specify station to process. **Note:** Only used if parameter *data* is wfdisc format.
    :type station: str
    :param byDay: Split data up by days into Stream (1 stream/day). This is useful if you want to avoid computer
                  memory issues. Note: Only used if parameter *data* is wfdisc format or directory.
    :type byDay: bool
    :param startdate: Specific startdate. Used to slice Stream for specific time range. Format ``YYYY-MM-DD HH:MM:SS``
    :type startdate: str
    :param enddate: Specific startdate. Used to slice Stream for specific time range. Format ``YYYY-MM-DD HH:MM:SS``
    :type enddate: str
    :param jul_start: Julian start date of data. This is only used when parameter *datatype* is ``"dir"`` or. It will
                      speed up file aggregtion, instead of looping and searching through a range of 365. It can also be
                      used for when parameter *datatype* is ``"wfdisc"`` for selecting which date to start running on
                      for a specific station.
    :type jul_start: int
    :param jul_end: Julian end date of data. This is only used when parameter *datatype* is ``"dir"``. It will speed up
                    file aggregtion, instead of looping and searching through a range of 365.
    :type jul_end: int
    :param generateMasks: Generate QC Masks. If True, masks will be generated when a value exceeds threshold
                          parameters.
    :type generateMasks: bool
    :param masksByTime: If masks are generated, the output mask will be in time if True, else it will be a binary mask.

    :type masksByTime: bool
    :param log: If filename is input, then a logger will be output, (``Default = "pycheron.log"``)
    :type log: str
    :param fortran: Use Fortran libs or not. **Note:** If libs will not compile or on a Windows Machine, set to False.
    :type fortran: bool
    :param database: If database filename defined, a db will be created (``Default = "pycheron.db"``)
    :type database: str
    :param session: Session name in database. Used for grouping different runs/experiments (``Default = "None"``).
                    **Note:** the defualt is a String ``"None"``, not ``None`` (python Null). Do **NOT** set it
                    ``None`` (python Null), otherwise it will not be able to filter correctly.
    :type session: str
    :param overwrite: If set to True, entries of the same snclq and start/end time will be overwritten.
    :type overwrite: bool
    :param to_csv: If set to true, database will output to csv files.
    :type to_csv: bool
    :param stationStartAt: **WFDISC ONLY** Specify the station to start processing at, instead of at the beginning of
                           the wfidisc file.
    :type stationStartAt: str

    **basicMetricStats**

    :param calcBasic: Calculates basicStatsMetric individually. **Note:** CalcAll must be set to ``False``.
    :type calcBasic: bool
    :param rmsThreshold: QC Mask threshold for RMS in basicStats.
    :type rmsThreshold: int
    :param maxThreshold: QC Mask threshold for Max in basicStats.
    :type maxThreshold: int
    :param minThreshold: QC Mask threshold for Min in basicStats.
    :type minThreshold: int
    :param medianThreshold: QC Mask threshold for Median in basicStats.
    :type medianThreshold: int
    :param meanThreshold: QC Mask threshold for Mean in basicStats.
    :type meanThreshold: int
    :param varianceThreshold: QC Mask threshold for Variance in basicStats.
    :type varianceThreshold: int
    :param stdThreshold: QC Mask threshold for STD in basicStats.
    :type stdThreshold: int

    **correlationMetric**

    :param calcCorr: Calculates correlationMetric individually. **Note:** CalcAll must be set to ``False``.
    :type calcCorr: bool

    **Parameters for crossCorrMetric**

    :param calcCrossCorr: Calculates crossCorrelationMetric individually. **Note:** CalcAll must be set to ``False``.

    :type calcCrossCorr: bool
    :param maxLagSecs: Maximum number of seconds of lag to use in crossCorrMetric. (Default = 10)
    :type maxLagSecs: int
    :param filt: Filter to use in crossCorrMetric. Current supported filters:

                * bandpass
                * bandstop
                * enevelope
                * highpass
                * integer_decimation
                * lowpass
                * lowpass_cheby_2

                For more in-depth explaination of filters. Please see :doc:`../pycheron.metrics.crossCorrMetric`
    :type filt: str
    :param freqmin: Depending on the filter specified in parameter *filt* for crossCorrMetric, either pass band low
                    corner frequency (bandpass, bandstop), filter corner frequency (lowpass), or the frequency
                    above which signals are attenuated with 95 dB.
    :type freqmin: int
    :param freqmax: Depending on the filter specified in :param filt for crossCorrMetric, either passband high
                    corner frequency (bandpass, bandstop), or filter corner frequency (highpass)
    :type freqmax: int
    :param corners: Filter corners/order for crossCorrMetric.
    :type corners: int
    :param zerophase: If True, applies filter once forward and once backwards in crossCorrMetric. This results
                      in twice the filter order but zero phase shift in the resulting filtered trace.

    :type zerophase: bool
    :param maxorder: Maximal order of designed cheby2 filter in crossCorrMetric. **Note:** option for lowpass_cheby2
                     only.
    :type maxorder: int
    :param ba: If True, returns only the filter coefficients (b,a) instead of filtering. **Note:** option for
               lowpass_cheby2 only.
    :type ba: bool
    :param freq_passband: If True, return additionally to the filtered data, then iteratively determine the
                          pass band frequency. **Note:** option for lowpass_cheby2 only.
    :type freq_passband: bool

    **DCOffSetTimesMetric**

    :param calcDcOffset: Calculates DCOffSetTimesMetric individually. **Note:** CalcAll must be set to ``False``.
    :type calcDcOffset: bool
    :param windowSecs: Chunk size (secs) used in DCOffset calculations
    :type windowSecs: int
    :param incrementSecs: Increment (secs) for starttime of sequential chunks (Default = windowSecs/2)
    :type incrementSecs: int
    :param threshold: Threshold used in the detection metric
    :type threshold: int

    **gapMetric**

    :param calcGap: Calculates gapMetric individually. **Note:** CalcAll must be set to ``False``.
    :type calcGap: bool
    :param separateMasks: If True, 2 masks will be create. 1 for gaps, and 1 for overlaps, if False these masks
                          will be merged and 3 masks will be returned. Gaps, Overlaps, Combined.
    :type separateMasks: bool
    :param completeDay: If True, if an individual day is does not start at 00:00:00 or end at 23:59:59, a small
                        trace will be made to complete the day. This should be used if, for example, have only
                        6 hours of data, but would like to consider that 6 hours out of a whole day, thus gaps
                        would be found and masks would be created around that segment of data.
    :type completeDay: bool

    **psdMetric**

    :param calcPsds: Calculates PSDS individually. **Note:** CalcAll must be set to ``False``.
    :type calcPsds: bool
    :param expLoPeriod: Low end of the period band used for calculating the exponential dead channel detector
                        (``Default= 4/trace sampling rate``)
    :type expLoPeriod: int
    :param expHiPeriod: High end of the period band used for calculating the exponential dead channel detector
    :type expHiPeriod: int
    :param linLoPeriod: Low end of the period band used for calculating the linear dead channel detector
                        (``Default= 16/trace sampling rate``)
    :type linLoPeriod: int
    :param linHiPeriod: High end of the period band used for calculating the linear dead channel detector
    :type linHiPeriod: int
    :param evalresp: Response information obtained from numpy.ndarray of freq, amp, phase information matching output
                     of ``client.evalresp``, optional for ``psdStatistics``. If ``None``, will pull from
                     ``client.evalresp``. Also accepts str path to RESP files
    :type evalresp: str or numpy.ndarray
    :param dcExpThreshold: Dead channel exponent threshold.
    :type dcExpThreshold: float
    :param pctBelowNoiseThreshold: Percent below NLNM threshold.
    :type pctBelowNoiseThreshold: int
    :param pctAboveNoiseThreshold: Percent above NLNM threshold.
    :type pctAboveNoiseThreshold: int
    :param dcLinThreshold: Dead channel linear threshold.
    :type dcLinThreshold: int
    :param num_gaps: Number of gaps threshold.
    :type num_gaps: int
    :param pctBelowNoiseThresholdRESP: Percent below NLNM for bad response threshold.
    :type pctBelowNoiseThresholdRESP: int
    :param pctAboveNoiseThresholdRESP: Percent above NLNM for bad response threshold.
    :type pctAboveNoiseThresholdRESP: int
    :param processesPSD: Number of processes to run. The max number of processes that can be run is
                         equal to the length of the Stream. Not recommended to go above 7 processes.
    :type processesPSD: int
    :param dcExpThresholdHour: dead channel exponent hourly threshold. Only used if ``byHourOn`` is ``True``
    :type dcExpThresholdHour: int
    :param dcLinThresholdHour: dead channel linear hourly threshold. Only used if ``byHouron`` set to ``True``
    :type dcLinThresholdHour: int
    :param byHourOn: Turn on hourly thresholding checks for ``dead_channel_exponential``, ``dead_channel_linear``,
                     and dead_channel_gsn checks
    :type byHourOn: bool

    **repeatedAmplitudeMetric**

    :param calcAmp: Calculates repeatedAmplitudeMetric individually. **Note:** CalcAll must be set to ``False``.
    :type calcAmp: bool
    :param minRep: minimum number of repeated adjacent values in a series, > 1.
    :type minRep: int

    **snrMetric**

    :param calcSNR: (bool) Calculates snrMetric individually. **Note:** CalcAll must be set to ``False``. (Default = False)
    :type calcSNR: bool
    :param algorithmSNR: A named algorithm to use for calculating SNR. Options include:

                        * 'splitWindow'
                        * 'staltaTrigger'
                        * 'pick'

                        Please see :doc:`../pycheron.metrics.snrMetric` for more detailed descriptions of algorithms.
    :type algorithmSNR: str
    :param windowSecsSNR: Width (secs) of the full window used in SNR calculation
    :type windowSecsSNR: int
    :param snrThreshold: SNR threshold at which if above, value would be masked.
    :type snrThreshold: int

    **sohMetric**

    :param calcSOH: Calculates sohMetric individually. **Note:** CalcAll must be set to ``False``.
    :type calcSOH: bool
    :param data_quality: If ``True``, only the data quality flags are returned.
    :type data_quality: bool
    :param activity: If ``True``, only the activity flags are returned.
    :type activity: bool
    :param io_clock: If ``True``, only the io/clock flags are returned.
    :type io_clock: bool

    **spikesMetric**

    :param calcSpikes: Calculates spikesMetric individually. **Note:** CalcAll must be set to ``False``.
    :type calcSpikes: bool
    :param windowSize: Size of window to be used (e.g., window size to roll over)
    :type windowSize: int
    :param thresholdSpikes: Initial threshold value for outlier detection
    :type thresholdSpikes: int
    :param selectivity: ``[0-1]`` Used in determining outliers, or ``None`` if ``fixedThreshold = True``
    :type selectivity: float
    :param fixedThreshold: If ``True``, set the threshold = threshold and ignore selectivity. If ``False``, then
                           the threshold is set to the maximum value of the hampel/hampelFilter function
                           output multiplied by the selectivity
    :type fixedThreshold: bool
    :param processesSpikes: Number of processes to bring up. Not recommended to go above 7 on personal machines.
    :type processesSpikes: int

    **staltaMetric**

    :param calcStalta: Calculates staltaMetric individually. **Note:** CalcAll must be set to ``False``.
    :type calcStalta: bool
    :param staSecs: Length of short term averaging window in secs
    :type staSecs: int
    :param ltaSecs: Length of long term averaging window in secs
    :type ltaSecs: int
    :param increment: Increment (in secs) used when sliding the average windows to the next location
    :type increment: int
    :param algorithmSTA: Algorithm to be used. Currently supported algorithms include (for more information see
                         references below and STALTA function docstrings):

                         * Classic_RR
                         * Classic_LR
                         * EarleAndShearer_envelope

    :type algorithmSTA: str

    **Parameters for deadChannelMetric**

    :param calcAllDeadChan: Calculates deadChannelMetric, deanChanADF, and deadChanMeanMetric individually.
                            **Note:** CalcAll must be set to ``False``.
    :type calcAllDeadChan: bool

    **deadChanADFMetric**

    :param dcADF_win_size: window size to increment the data by in samples, each window will have a ADF run on it
                           (``DEFAULT = 500``, assumes 100sps). If None, window size determined by sample rate of each
                           trace in the stream object
    :type dcADF_win_size: int
    :param dcADF_pval_thresh: p value threshold to trigger ADF rejection of null hypothesis (e.g., time series is
                              assumed stationary with rejection of null hypothesis)
    :type dcADF_pval_thresh: float
    :param dcADF_threshold: threshold to trigger ADF rejection of null hypothesis using difference between DF test statistic
                            and 1% critical value
    :type dcADF_threshold: float
    :param dcADF_use_thresh: String option for choosing whether to use the pvalue_threshold or the threshold for
                             determining whether null hypothesis fails (DEFAULT = pvalue). Options are pvalue or
                             diffs.
    :type dcADF_use_thresh: float

    **deadChanMeanMetric**

    :param dcMean_win_size: number of samples to use for windowing. When ``None``, determines the number
                            of samples to use based on the sample rate
    :type dcMean_win_size: int
    :param dcMean_thresh: threshold to trigger deadChan for within rolling window. Recommended
                          to keep this default value for triggered a dead channel
    :type dcMean_thresh: float

    **transferFunctionMetric**

    :param calcTransfer: Calculates transferFunctionMetric individually. **Note:** CalcAll must be set to ``False``.
    :type calcTransfer: bool

    **calibrationMetric**

    :param calcCal: Calculates scalibrationMetric individually. **Note:** CalcAll must be set to ``False``.
    :type calcCal: bool
    :param cal_metric_store: (Default = None)

    **plots**

    :param plots: (bool) If True, plots will be created. This will slow down the program considerably. It is
                         recommended to use the bulk plotting functions once calculations are complete.
    :type plots: bool
    :param pdfModel: Noise Model. Supports:

                     * ``"gsn"`` - Global Seismic Network 2004
                     * ``"nlnm"`` - New Low Noise Model

    :type pdfModel: str
    :param per_arr: Period array for PDF Mode Colorline plot. If none, periods will be:
                    ``[.1, 1.0, 6.5, 10, 30.8, 103.7]``.
    :type per_arr: numpy.array or list
    :param showNoiseModel: If ``True``, display noise models on the plot
    :type showNoiseModel: bool
    :param showMaxMin: If ``True``, display max and min on the plot
    :type showMaxMin: bool
    :param showMode: If ``True``, display mode on the plot
    :type showMode: bool
    :param showMean: If ``True``, display mean on the plot
    :type showMean: bool
    :param showMedian: If ``True``, display median on the plot
    :type showMedian: bool
    :param showEnvelope: If ``True``, percentile envelope will be plotted
    :type showEnvelope: bool
    :param envelopeType: Envelope to be plotted. Options:

                               * ``"10_90"`` - 10th and 90th percentile envelope
                               * ``"05_95"`` - 5th annd 95th percentile envelope

    :type envelopeType: bool
    :param showSingle: If ``True``, plot single percentile.
    :type showSingle: bool
    :param singleType: Single line percentile to be plotted. Options:

                             * ``"5"`` - 5th percentile
                             * ``"10"`` - 10th percentile
                             * ``"90"`` - 90th percentile
                             * ``"95"`` - 95th percentile

    :type singleType: bool
    :param ylo: min dB to show on the y-axis scale
    :type ylo: int
    :param yhi: max dB to show on the y-axis scale
    :type yhi: int
    :param min_stations: Minimum number of stations to have calculated before Station Ranking Plots will be
                         plotted.
    :type min_stations: int
    :param rank_by: Array of periods to rank by in Station Ranking Plots. If None, periods will be
                    ``[1, 6.5, 30, 100]``
    :type rank_by: numpy.array or list
    :param timespan: Select the timespan for psd/pdf plots. Options are ``1`` (week), ``2`` (weeks), ``4`` (weeks).
                     Any other option will defualt to plotting all psds (``Default="All"``)
    :type timespan: int or string
    """

    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    if database is not None:
        database = output_dir + "/" + database
        db = Database(database, session, overwrite)

    # creating logger
    if log is not None:
        logger = Logger(log)

    start = time.time()
    # for single .mseed file
    if datatype == "mseed":
        # read in file
        st_raw = obspy.read(data, format="MSEED")
        # set julday to "All", means processing all data, not by individual days
        julday = "All"
        # getting all unique stations from Stream
        un_sta = getUniqSta(st_raw)
        # looping through station list and calling metrics on each station
        for i in un_sta:
            st_sta = getSta(st_raw, i)
            _call_pycheron_wrapper(
                output_dir,
                st_sta,
                calcAll,
                calcPsds,
                calcBasic,
                calcCorr,
                calcCrossCorr,
                calcGap,
                calcAmp,
                calcSNR,
                calcSOH,
                calcStalta,
                calcDcOffset,
                calcTransfer,
                calcSpikes,
                calcAllDeadChan,
                calcCal,
                julday,
                startdate,
                enddate,
                jul_start,
                jul_end,
                generateMasks,
                masksByTime,
                rmsThreshold,
                maxThreshold,
                minThreshold,
                medianThreshold,
                meanThreshold,
                varianceThreshold,
                stdThreshold,
                maxLagSecs,
                filt,
                freqmin,
                freqmax,
                corners,
                zerophase,
                maxorder,
                ba,
                freq_passband,
                windowSecs,
                incrementSecs,
                threshold,
                separateMasks,
                expLoPeriod,
                expHiPeriod,
                linLoPeriod,
                linHiPeriod,
                evalresp,
                dcExpThreshold,
                pctBelowNoiseThreshold,
                pctAboveNoiseThreshold,
                dcLinThreshold,
                num_gaps,
                pctBelowNoiseThresholdRESP,
                pctAboveNoiseThresholdRESP,
                minRep,
                algorithmSNR,
                windowSecsSNR,
                snrThreshold,
                data_quality,
                activity,
                io_clock,
                windowSize,
                thresholdSpikes,
                selectivity,
                fixedThreshold,
                staSecs,
                ltaSecs,
                increment,
                algorithmSTA,
                plots,
                pdfModel,
                per_arr,
                showNoiseModel,
                showMaxMin,
                showMode,
                showMean,
                showMedian,
                showEnvelope,
                envelopeType,
                showSingle,
                singleType,
                ylo,
                yhi,
                min_stations,
                rank_by,
                processesPSD,
                processesSpikes,
                logger,
                fortran,
                timespan,
                dcADF_win_size,
                dcADF_threshold,
                dcADF_pval_thresh,
                dcADF_use_thresh,
                dcMean_win_size,
                dcMean_thresh,
                cal_metric_store,
                dcExpThresholdHour,
                dcLinThresholdHour,
                byHourOn,
                db,
                session,
            )

    # for Obspy Stream object -> mostly for SAPL integration
    elif datatype == "stream":
        st_raw = data
        # set julday to "All", means processing all data, not by individual days
        julday = "All"
        # getting all unique stations from Stream
        un_sta = getUniqSta(st_raw)
        # looping through station list and calling metrics on each station
        for i in un_sta:
            st_sta = getSta(st_raw, i)
            _call_pycheron_wrapper(
                output_dir,
                st_sta,
                calcAll,
                calcPsds,
                calcBasic,
                calcCorr,
                calcCrossCorr,
                calcGap,
                calcAmp,
                calcSNR,
                calcSOH,
                calcStalta,
                calcDcOffset,
                calcTransfer,
                calcSpikes,
                calcAllDeadChan,
                calcCal,
                julday,
                startdate,
                enddate,
                jul_start,
                jul_end,
                generateMasks,
                masksByTime,
                rmsThreshold,
                maxThreshold,
                minThreshold,
                medianThreshold,
                meanThreshold,
                varianceThreshold,
                stdThreshold,
                maxLagSecs,
                filt,
                freqmin,
                freqmax,
                corners,
                zerophase,
                maxorder,
                ba,
                freq_passband,
                windowSecs,
                incrementSecs,
                threshold,
                separateMasks,
                expLoPeriod,
                expHiPeriod,
                linLoPeriod,
                linHiPeriod,
                evalresp,
                dcExpThreshold,
                pctBelowNoiseThreshold,
                pctAboveNoiseThreshold,
                dcLinThreshold,
                num_gaps,
                pctBelowNoiseThresholdRESP,
                pctAboveNoiseThresholdRESP,
                minRep,
                algorithmSNR,
                windowSecsSNR,
                snrThreshold,
                data_quality,
                activity,
                io_clock,
                windowSize,
                thresholdSpikes,
                selectivity,
                fixedThreshold,
                staSecs,
                ltaSecs,
                increment,
                algorithmSTA,
                plots,
                pdfModel,
                per_arr,
                showNoiseModel,
                showMaxMin,
                showMode,
                showMean,
                showMedian,
                showEnvelope,
                envelopeType,
                showSingle,
                singleType,
                ylo,
                yhi,
                min_stations,
                rank_by,
                processesPSD,
                processesSpikes,
                logger,
                fortran,
                timespan,
                dcADF_win_size,
                dcADF_threshold,
                dcADF_pval_thresh,
                dcADF_use_thresh,
                dcMean_win_size,
                dcMean_thresh,
                cal_metric_store,
                dcExpThresholdHour,
                dcLinThresholdHour,
                byHourOn,
                db,
                session,
            )

    # for CSS/S4 Wfdisc file
    elif datatype == "wfdisc":
        # gets all unique stations from wfdisc file
        un_sta = get_wfdisc_stations(data)

        # if station not specified, process all stations in file
        if station is None:
            if stationStartAt is not None:
                idx = list(un_sta).index(stationStartAt)
                un_sta = un_sta[idx:]
            for i in un_sta:
                # convert waveform files into Obspy Stream objects
                st_sta = css2stream(data, network, i, byDay)
                # if stream to be entire station
                if byDay == False:
                    julday = "All"
                    _call_pycheron_wrapper(
                        output_dir,
                        st_sta,
                        calcAll,
                        calcPsds,
                        calcBasic,
                        calcCorr,
                        calcCrossCorr,
                        calcGap,
                        calcAmp,
                        calcSNR,
                        calcSOH,
                        calcStalta,
                        calcDcOffset,
                        calcTransfer,
                        calcSpikes,
                        calcAllDeadChan,
                        calcCal,
                        julday,
                        startdate,
                        enddate,
                        jul_start,
                        jul_end,
                        generateMasks,
                        masksByTime,
                        rmsThreshold,
                        maxThreshold,
                        minThreshold,
                        medianThreshold,
                        meanThreshold,
                        varianceThreshold,
                        stdThreshold,
                        maxLagSecs,
                        filt,
                        freqmin,
                        freqmax,
                        corners,
                        zerophase,
                        maxorder,
                        ba,
                        freq_passband,
                        windowSecs,
                        incrementSecs,
                        threshold,
                        separateMasks,
                        expLoPeriod,
                        expHiPeriod,
                        linLoPeriod,
                        linHiPeriod,
                        evalresp,
                        dcExpThreshold,
                        pctBelowNoiseThreshold,
                        pctAboveNoiseThreshold,
                        dcLinThreshold,
                        num_gaps,
                        pctBelowNoiseThresholdRESP,
                        pctAboveNoiseThresholdRESP,
                        minRep,
                        algorithmSNR,
                        windowSecsSNR,
                        snrThreshold,
                        data_quality,
                        activity,
                        io_clock,
                        windowSize,
                        thresholdSpikes,
                        selectivity,
                        fixedThreshold,
                        staSecs,
                        ltaSecs,
                        increment,
                        algorithmSTA,
                        plots,
                        pdfModel,
                        per_arr,
                        showNoiseModel,
                        showMaxMin,
                        showMode,
                        showMean,
                        showMedian,
                        showEnvelope,
                        envelopeType,
                        showSingle,
                        singleType,
                        ylo,
                        yhi,
                        min_stations,
                        rank_by,
                        processesPSD,
                        processesSpikes,
                        logger,
                        fortran,
                        timespan,
                        dcADF_win_size,
                        dcADF_threshold,
                        dcADF_pval_thresh,
                        dcADF_use_thresh,
                        dcMean_win_size,
                        dcMean_thresh,
                        cal_metric_store,
                        dcExpThresholdHour,
                        dcLinThresholdHour,
                        byHourOn,
                        db,
                        session,
                    )
                else:
                    # if stream to be split up by days, loop through each day in stream
                    for st in st_sta:
                        if (
                            jul_start is None
                            or st[0].stats.starttime.julday >= jul_start
                        ):
                            julday = st[0].stats.starttime.julday
                            _call_pycheron_wrapper(
                                output_dir,
                                st,
                                calcAll,
                                calcPsds,
                                calcBasic,
                                calcCorr,
                                calcCrossCorr,
                                calcGap,
                                calcAmp,
                                calcSNR,
                                calcSOH,
                                calcStalta,
                                calcDcOffset,
                                calcTransfer,
                                calcSpikes,
                                calcAllDeadChan,
                                calcCal,
                                julday,
                                startdate,
                                enddate,
                                jul_start,
                                jul_end,
                                generateMasks,
                                masksByTime,
                                rmsThreshold,
                                maxThreshold,
                                minThreshold,
                                medianThreshold,
                                meanThreshold,
                                varianceThreshold,
                                stdThreshold,
                                maxLagSecs,
                                filt,
                                freqmin,
                                freqmax,
                                corners,
                                zerophase,
                                maxorder,
                                ba,
                                freq_passband,
                                windowSecs,
                                incrementSecs,
                                threshold,
                                separateMasks,
                                expLoPeriod,
                                expHiPeriod,
                                linLoPeriod,
                                linHiPeriod,
                                evalresp,
                                dcExpThreshold,
                                pctBelowNoiseThreshold,
                                pctAboveNoiseThreshold,
                                dcLinThreshold,
                                num_gaps,
                                pctBelowNoiseThresholdRESP,
                                pctAboveNoiseThresholdRESP,
                                minRep,
                                algorithmSNR,
                                windowSecsSNR,
                                snrThreshold,
                                data_quality,
                                activity,
                                io_clock,
                                windowSize,
                                thresholdSpikes,
                                selectivity,
                                fixedThreshold,
                                staSecs,
                                ltaSecs,
                                increment,
                                algorithmSTA,
                                plots,
                                pdfModel,
                                per_arr,
                                showNoiseModel,
                                showMaxMin,
                                showMode,
                                showMean,
                                showMedian,
                                showEnvelope,
                                envelopeType,
                                showSingle,
                                singleType,
                                ylo,
                                yhi,
                                min_stations,
                                rank_by,
                                processesPSD,
                                processesSpikes,
                                logger,
                                fortran,
                                timespan,
                                dcADF_win_size,
                                dcADF_threshold,
                                dcADF_pval_thresh,
                                dcADF_use_thresh,
                                dcMean_win_size,
                                dcMean_thresh,
                                cal_metric_store,
                                dcExpThresholdHour,
                                dcLinThresholdHour,
                                byHourOn,
                                db,
                                session,
                            )
        # if specific Station is specified
        else:
            # convert waveform file to Stream
            st_sta = css2stream(data, network, station, byDay)
            # if stream to be entire station
            if byDay == False:
                julday = "All"
                _call_pycheron_wrapper(
                    output_dir,
                    st_sta,
                    calcAll,
                    calcPsds,
                    calcBasic,
                    calcCorr,
                    calcCrossCorr,
                    calcGap,
                    calcAmp,
                    calcSNR,
                    calcSOH,
                    calcStalta,
                    calcDcOffset,
                    calcTransfer,
                    calcSpikes,
                    calcAllDeadChan,
                    calcCal,
                    julday,
                    startdate,
                    enddate,
                    jul_start,
                    jul_end,
                    generateMasks,
                    masksByTime,
                    rmsThreshold,
                    maxThreshold,
                    minThreshold,
                    medianThreshold,
                    meanThreshold,
                    varianceThreshold,
                    stdThreshold,
                    maxLagSecs,
                    filt,
                    freqmin,
                    freqmax,
                    corners,
                    zerophase,
                    maxorder,
                    ba,
                    freq_passband,
                    windowSecs,
                    incrementSecs,
                    threshold,
                    separateMasks,
                    expLoPeriod,
                    expHiPeriod,
                    linLoPeriod,
                    linHiPeriod,
                    evalresp,
                    dcExpThreshold,
                    pctBelowNoiseThreshold,
                    pctAboveNoiseThreshold,
                    dcLinThreshold,
                    num_gaps,
                    pctBelowNoiseThresholdRESP,
                    pctAboveNoiseThresholdRESP,
                    minRep,
                    algorithmSNR,
                    windowSecsSNR,
                    snrThreshold,
                    data_quality,
                    activity,
                    io_clock,
                    windowSize,
                    thresholdSpikes,
                    selectivity,
                    fixedThreshold,
                    staSecs,
                    ltaSecs,
                    increment,
                    algorithmSTA,
                    plots,
                    pdfModel,
                    per_arr,
                    showNoiseModel,
                    showMaxMin,
                    showMode,
                    showMean,
                    showMedian,
                    showEnvelope,
                    envelopeType,
                    showSingle,
                    singleType,
                    ylo,
                    yhi,
                    min_stations,
                    rank_by,
                    processesPSD,
                    processesSpikes,
                    logger,
                    fortran,
                    timespan,
                    dcADF_win_size,
                    dcADF_threshold,
                    dcADF_pval_thresh,
                    dcADF_use_thresh,
                    dcMean_win_size,
                    dcMean_thresh,
                    cal_metric_store,
                    dcExpThresholdHour,
                    dcLinThresholdHour,
                    byHourOn,
                    db,
                    session,
                )
            # if stream to be split up by days, loop through each day in stream
            else:
                for st in st_sta:
                    if jul_start is None or st[0].stats.starttime.julday >= jul_start:
                        julday = st[0].stats.starttime.julday
                        _call_pycheron_wrapper(
                            output_dir,
                            st,
                            calcAll,
                            calcPsds,
                            calcBasic,
                            calcCorr,
                            calcCrossCorr,
                            calcGap,
                            calcAmp,
                            calcSNR,
                            calcSOH,
                            calcStalta,
                            calcDcOffset,
                            calcTransfer,
                            calcSpikes,
                            calcAllDeadChan,
                            calcCal,
                            julday,
                            startdate,
                            enddate,
                            jul_start,
                            jul_end,
                            generateMasks,
                            masksByTime,
                            rmsThreshold,
                            maxThreshold,
                            minThreshold,
                            medianThreshold,
                            meanThreshold,
                            varianceThreshold,
                            stdThreshold,
                            maxLagSecs,
                            filt,
                            freqmin,
                            freqmax,
                            corners,
                            zerophase,
                            maxorder,
                            ba,
                            freq_passband,
                            windowSecs,
                            incrementSecs,
                            threshold,
                            separateMasks,
                            expLoPeriod,
                            expHiPeriod,
                            linLoPeriod,
                            linHiPeriod,
                            evalresp,
                            dcExpThreshold,
                            pctBelowNoiseThreshold,
                            pctAboveNoiseThreshold,
                            dcLinThreshold,
                            num_gaps,
                            pctBelowNoiseThresholdRESP,
                            pctAboveNoiseThresholdRESP,
                            minRep,
                            algorithmSNR,
                            windowSecsSNR,
                            snrThreshold,
                            data_quality,
                            activity,
                            io_clock,
                            windowSize,
                            thresholdSpikes,
                            selectivity,
                            fixedThreshold,
                            staSecs,
                            ltaSecs,
                            increment,
                            algorithmSTA,
                            plots,
                            pdfModel,
                            per_arr,
                            showNoiseModel,
                            showMaxMin,
                            showMode,
                            showMean,
                            showMedian,
                            showEnvelope,
                            envelopeType,
                            showSingle,
                            singleType,
                            ylo,
                            yhi,
                            min_stations,
                            rank_by,
                            processesPSD,
                            processesSpikes,
                            logger,
                            fortran,
                            timespan,
                            dcADF_win_size,
                            dcADF_threshold,
                            dcADF_pval_thresh,
                            dcADF_use_thresh,
                            dcMean_win_size,
                            dcMean_thresh,
                            cal_metric_store,
                            dcExpThresholdHour,
                            dcLinThresholdHour,
                            byHourOn,
                            db,
                            session,
                        )

    elif datatype == "dir":
        # create array of jul days
        if jul_start and jul_end is not None:
            date_range = np.arange(jul_start, jul_end + 1)
        else:
            # if jul date range is None, then process all of the files
            date_range = np.arange(1, 366, 1)

        # Convert julian days to real julian day numbers for those less than 100. The reason we do this is that
        # np.arange will not create an appropriate date_range if inputs are 001, 032 for example. So we take the user
        # input and convert to proper julian days so that we can ensure that the correct data files are being grabbed
        new_dates = []
        for i in date_range:
            if i < 100:
                if i < 10:
                    i = str(i)
                    i = "00" + i
                    new_dates.append(i)
                elif i >= 10 and i < 100:
                    i = str(i)
                    i = "0" + i
                    new_dates.append(i)
            else:
                new_dates.append(str(i))

        # Create empty sta_list and chan_list for looping
        sta_list = []
        chan_list = []
        # looping through files
        for i in os.listdir(data):
            # if file ends with jul date add to list
            if i.endswith(".mseed") or i.endswith(".w"):
                # get station from name
                name = i.split("_")[0]
                sta_list.append(name)
                chan_list.append(i.split("_")[1])
        # get unique stations
        un_chans = set(chan_list)
        un_chans = list(un_chans)
        un_chans.sort()
        un_sta = set(sta_list)
        un_sta = list(un_sta)
        un_sta.sort()

        # aggregate all days of one station into single stream
        if byDay == False:
            # looping through files
            for i in un_sta:
                stream = []
                for j in new_dates:
                    for k in os.listdir(data):
                        # if i matches station
                        if fnmatch.fnmatch(k, "*" + i + "*") and k.endswith(
                            j + ".mseed"
                        ):
                            # st filename of specific station
                            st = data + "/" + str(k)
                            stream.append(st)
                stream = set(stream)

                # initialize first stream
                st_sta = obspy.read(list(stream)[0])
                if len(stream) > 1:
                    for s in range(1, len(stream)):
                        st_sta = st_sta + obspy.read(list(stream)[s])
                print("Finished Aggregating Streams for: " + i)
                julday = "All"
                _call_pycheron_wrapper(
                    output_dir,
                    st_sta,
                    calcAll,
                    calcPsds,
                    calcBasic,
                    calcCorr,
                    calcCrossCorr,
                    calcGap,
                    calcAmp,
                    calcSNR,
                    calcSOH,
                    calcStalta,
                    calcDcOffset,
                    calcTransfer,
                    calcSpikes,
                    calcAllDeadChan,
                    calcCal,
                    julday,
                    startdate,
                    enddate,
                    jul_start,
                    jul_end,
                    generateMasks,
                    masksByTime,
                    rmsThreshold,
                    maxThreshold,
                    minThreshold,
                    medianThreshold,
                    meanThreshold,
                    varianceThreshold,
                    stdThreshold,
                    maxLagSecs,
                    filt,
                    freqmin,
                    freqmax,
                    corners,
                    zerophase,
                    maxorder,
                    ba,
                    freq_passband,
                    windowSecs,
                    incrementSecs,
                    threshold,
                    separateMasks,
                    expLoPeriod,
                    expHiPeriod,
                    linLoPeriod,
                    linHiPeriod,
                    evalresp,
                    dcExpThreshold,
                    pctBelowNoiseThreshold,
                    pctAboveNoiseThreshold,
                    dcLinThreshold,
                    num_gaps,
                    pctBelowNoiseThresholdRESP,
                    pctAboveNoiseThresholdRESP,
                    minRep,
                    algorithmSNR,
                    windowSecsSNR,
                    snrThreshold,
                    data_quality,
                    activity,
                    io_clock,
                    windowSize,
                    thresholdSpikes,
                    selectivity,
                    fixedThreshold,
                    staSecs,
                    ltaSecs,
                    increment,
                    algorithmSTA,
                    plots,
                    pdfModel,
                    per_arr,
                    showNoiseModel,
                    showMaxMin,
                    showMode,
                    showMean,
                    showMedian,
                    showEnvelope,
                    envelopeType,
                    showSingle,
                    singleType,
                    ylo,
                    yhi,
                    min_stations,
                    rank_by,
                    processesPSD,
                    processesSpikes,
                    logger,
                    fortran,
                    timespan,
                    dcADF_win_size,
                    dcADF_threshold,
                    dcADF_pval_thresh,
                    dcADF_use_thresh,
                    dcMean_win_size,
                    dcMean_thresh,
                    cal_metric_store,
                    dcExpThresholdHour,
                    dcLinThresholdHour,
                    byHourOn,
                    db,
                    session,
                )
        # process by individual day
        else:
            # loop through station
            for i in un_sta:
                stream = []
                date_list = []
                for j in new_dates:
                    day = []
                    for k in os.listdir(data):
                        # if i matches station
                        if fnmatch.fnmatch(k, "*" + i + "*") and k.endswith(
                            j + ".mseed"
                        ):
                            # st filename of specific station
                            st = data + "/" + str(k)
                            day.append(st)
                            date_list.append(j)
                    day = set(day)
                    date = set(date_list)
                    if day:
                        stream.append(day)
                date = sorted(date)
                for k in range(len(stream)):
                    # initialize first stream
                    julday = date[k]
                    files = stream[k]
                    st_sta = obspy.read(list(files)[0])
                    if len(files) > 1:
                        for s in range(1, len(files)):
                            st_sta = st_sta + obspy.read(list(files)[s])
                    _call_pycheron_wrapper(
                        output_dir,
                        st_sta,
                        calcAll,
                        calcPsds,
                        calcBasic,
                        calcCorr,
                        calcCrossCorr,
                        calcGap,
                        calcAmp,
                        calcSNR,
                        calcSOH,
                        calcStalta,
                        calcDcOffset,
                        calcTransfer,
                        calcSpikes,
                        calcAllDeadChan,
                        calcCal,
                        julday,
                        startdate,
                        enddate,
                        jul_start,
                        jul_end,
                        generateMasks,
                        masksByTime,
                        rmsThreshold,
                        maxThreshold,
                        minThreshold,
                        medianThreshold,
                        meanThreshold,
                        varianceThreshold,
                        stdThreshold,
                        maxLagSecs,
                        filt,
                        freqmin,
                        freqmax,
                        corners,
                        zerophase,
                        maxorder,
                        ba,
                        freq_passband,
                        windowSecs,
                        incrementSecs,
                        threshold,
                        separateMasks,
                        expLoPeriod,
                        expHiPeriod,
                        linLoPeriod,
                        linHiPeriod,
                        evalresp,
                        dcExpThreshold,
                        pctBelowNoiseThreshold,
                        pctAboveNoiseThreshold,
                        dcLinThreshold,
                        num_gaps,
                        pctBelowNoiseThresholdRESP,
                        pctAboveNoiseThresholdRESP,
                        minRep,
                        algorithmSNR,
                        windowSecsSNR,
                        snrThreshold,
                        data_quality,
                        activity,
                        io_clock,
                        windowSize,
                        thresholdSpikes,
                        selectivity,
                        fixedThreshold,
                        staSecs,
                        ltaSecs,
                        increment,
                        algorithmSTA,
                        plots,
                        pdfModel,
                        per_arr,
                        showNoiseModel,
                        showMaxMin,
                        showMode,
                        showMean,
                        showMedian,
                        showEnvelope,
                        envelopeType,
                        showSingle,
                        singleType,
                        ylo,
                        yhi,
                        min_stations,
                        rank_by,
                        processesPSD,
                        processesSpikes,
                        logger,
                        fortran,
                        timespan,
                        dcADF_win_size,
                        dcADF_threshold,
                        dcADF_pval_thresh,
                        dcADF_use_thresh,
                        dcMean_win_size,
                        dcMean_thresh,
                        cal_metric_store,
                        dcExpThresholdHour,
                        dcLinThresholdHour,
                        byHourOn,
                        db,
                        session,
                    )
        timestart = (time.time() - start) / 60
        logger.log(
            "callPycheron(): Total program runtime in minutes - " + str(timestart)
        )

    else:
        logger.error(
            "callPycheron(): Must provide either .mseed, .wfdisc, stream object or data directory"
        )
        sys.exit(
            "callPycheron(): Must provide either .mseed, .wfdisc, stream object or data directory"
        )

    if to_csv:
        db = Database(database)
        db.to_csv(output_dir + "/pycheron")


# -----------------------------------------------------------------------------
# Pycheron Wrapper Code
# -----------------------------------------------------------------------------


def _call_pycheron_wrapper(
    output_dir,
    st_sta,
    calcAll,
    calcPsds,
    calcBasic,
    calcCorr,
    calcCrossCorr,
    calcGap,
    calcAmp,
    calcSNR,
    calcSOH,
    calcStalta,
    calcDcOffset,
    calcTransfer,
    calcSpikes,
    calcAllDeadChan,
    calcCal,
    julday,
    startdate,
    enddate,
    jul_start,
    jul_end,
    generateMasks,
    masksByTime,
    rmsThreshold,
    maxThreshold,
    minThreshold,
    medianThreshold,
    meanThreshold,
    varianceThreshold,
    stdThreshold,
    maxLagSecs,
    filt,
    freqmin,
    freqmax,
    corners,
    zerophase,
    maxorder,
    ba,
    freq_passband,
    windowSecs,
    incrementSecs,
    threshold,
    separateMasks,
    expLoPeriod,
    expHiPeriod,
    linLoPeriod,
    linHiPeriod,
    evalresp,
    dcExpThreshold,
    pctBelowNoiseThreshold,
    pctAboveNoiseThreshold,
    dcLinThreshold,
    num_gaps,
    pctBelowNoiseThresholdRESP,
    pctAboveNoiseThresholdRESP,
    minRep,
    algorithmSNR,
    windowSecsSNR,
    snrThreshold,
    data_quality,
    activity,
    io_clock,
    windowSize,
    thresholdSpikes,
    selectivity,
    fixedThreshold,
    staSecs,
    ltaSecs,
    increment,
    algorithmSTA,
    plots,
    pdfModel,
    per_arr,
    showNoiseModel,
    showMaxMin,
    showMode,
    showMean,
    showMedian,
    showEnvelope,
    envelopeType,
    showSingle,
    singleType,
    ylo,
    yhi,
    min_stations,
    rank_by,
    processesPSD,
    processesSpikes,
    logger,
    fortran,
    timespan,
    dcADF_win_size,
    dcADF_threshold,
    dcADF_pval_thresh,
    dcADF_use_thresh,
    dcMean_win_size,
    dcMean_thresh,
    cal_metric_store,
    dcExpThresholdHour,
    dcLinThresholdHour,
    byHourOn,
    database,
    session,
):
    """
    Internal function that callPycheron calls. This does the majority of the processing.
    """
    st = st_sta.copy()
    date = str(julday)
    # split trace into time range
    if startdate and enddate is not None:
        starttime = UTCDateTime(startdate)
        endtime = UTCDateTime(enddate)
        st = st.slice(starttime, endtime)

    # merges traces needed for crossCorr and corrMetric
    stM = st_sta.copy()
    stM = stM.merge()

    # getting snclq for network and station
    tr = st[0]
    snclq = tr.get_id()
    snclq = snclq.split(".")

    # network dir
    network = snclq[0]
    net_dir = output_dir + "/" + network
    if os.path.exists(net_dir) == False:
        os.mkdir(net_dir)

    # station dir
    station = snclq[1]
    sta_dir = net_dir + "/" + station
    if os.path.exists(sta_dir) == False:
        os.mkdir(sta_dir)

    # channel dirs
    chans = []
    for i in range(len(st)):
        tr = st[i]
        snclq = tr.get_id()
        channel = getChannelName(snclq)
        chans.append(channel)

    # get unique channels from stream
    un_chans = np.unique(chans)

    print("-----------------------------------------------")
    print(
        (
            "Beginning Metric Calculations for "
            + network
            + "."
            + station
            + " Jul Date: "
            + date
        )
    )
    print("-----------------------------------------------")
    # creating PSD folder for Station

    # # ----------------------psdMetric---------------------------------------------
    if calcAll or calcPsds:
        st_psd = st.copy()

        _psd_wrapper(
            st_psd,
            expLoPeriod,
            expHiPeriod,
            linLoPeriod,
            linHiPeriod,
            evalresp,
            generateMasks,
            masksByTime,
            dcExpThreshold,
            pctBelowNoiseThreshold,
            pctAboveNoiseThreshold,
            rmsThreshold,
            dcLinThreshold,
            num_gaps,
            pctBelowNoiseThresholdRESP,
            pctAboveNoiseThresholdRESP,
            processesPSD,
            network,
            station,
            logger,
            byHourOn,
            database,
        )

    # # ----------------------basicStatsMetric--------------------------------------
    if calcAll or calcBasic:
        st_basic = st.copy()
        t_basic = Process(
            target=_basic_stats_wrapper,
            args=(
                st_basic,
                rmsThreshold,
                maxThreshold,
                minThreshold,
                medianThreshold,
                meanThreshold,
                varianceThreshold,
                stdThreshold,
                generateMasks,
                masksByTime,
                network,
                station,
                logger,
                database,
            ),
        )
        t_basic.start()

    # # ----------------------correlationMetric-------------------------------------
    if calcAll or calcCorr:
        tmp = 0  # needed otherwise Process will think st object a mappable list.
        t_corr = Process(
            target=_corr_wrapper, args=(stM, network, station, logger, database)
        )
        t_corr.start()

    # # ----------------------crossCorrMetric---------------------------------------
    if calcAll or calcCrossCorr:
        st_cross = stM.copy()
        t_cross = Process(
            target=_cross_corr_wrapper,
            args=(
                stM,
                freqmin,
                maxLagSecs,
                filt,
                freqmax,
                corners,
                zerophase,
                maxorder,
                ba,
                freq_passband,
                network,
                station,
                logger,
                database,
            ),
        )
        t_cross.start()

    # # ----------------------DCOffSetTimesMetric-----------------------------------
    if calcAll or calcDcOffset:
        st_offset = st.copy()
        t_offset = Process(
            target=_offset_wrapper,
            args=(
                st_offset,
                generateMasks,
                masksByTime,
                windowSecs,
                incrementSecs,
                threshold,
                network,
                station,
                logger,
                database,
            ),
        )
        t_offset.start()

    # # ----------------------gapMetric---------------------------------------------
    if calcAll or calcGap:
        st_gap = st.copy()
        t_gap = Process(
            target=_gap_wrapper,
            args=(
                st_gap,
                generateMasks,
                separateMasks,
                masksByTime,
                network,
                station,
                logger,
                database,
            ),
        )  # calculating metric
        t_gap.start()

    # # ----------------------repeatedAmplitude-------------------------------------
    if calcAll or calcAmp:
        st_amp = st.copy()
        t_amp = Process(
            target=_rep_amp_wrapper,
            args=(
                st_amp,
                minRep,
                generateMasks,
                network,
                station,
                logger,
                False,
                database,
            ),
        )
        t_amp.start()

    # # ----------------------snrMetric---------------------------------------------
    if calcAll or calcSNR:
        st_snr = st.copy()
        t_snr = Process(
            target=_snr_wrapper,
            args=(
                st_snr,
                algorithmSNR,
                windowSecsSNR,
                snrThreshold,
                generateMasks,
                masksByTime,
                network,
                station,
                logger,
                database,
            ),
        )
        t_snr.start()

    # # ----------------------sohMetric---------------------------------------------
    if calcAll or calcSOH:
        # creating temp MSEED file because Obspy function used within sohMetric, requires actual file, not St object

        t_soh = Process(
            target=_soh_wrapper,
            args=(
                st,
                data_quality,
                activity,
                io_clock,
                network,
                station,
                logger,
                database,
            ),
        )
        t_soh.start()

    # # ----------------------staltaMetric------------------------------------------
    if calcAll or calcStalta:
        st_stalta = st.copy()
        t_stalta = Process(
            target=_stalta_wrapper,
            args=(
                st_stalta,
                staSecs,
                ltaSecs,
                increment,
                algorithmSTA,
                network,
                station,
                logger,
                fortran,
                database,
            ),
        )
        t_stalta.start()

    # # ----------------------transferFunctionMetric--------------------------------
    if calcAll or calcTransfer:
        tmp = 0  # needed otherwise Process will think st object a mappable list.
        st_trans = st.copy()
        t_transfer = Process(
            target=_transfer_wrapper,
            args=(st_trans, network, station, logger, database),
        )
        t_transfer.start()

    # # ----------------------spikesMetric------------------------------------------
    if calcAll or calcSpikes:
        st_spikes = st.copy()
        t_spikes = Process(
            target=_spikes_wrapper,
            args=(
                st_spikes,
                windowSize,
                thresholdSpikes,
                selectivity,
                fixedThreshold,
                processesSpikes,
                network,
                station,
                logger,
                fortran,
                generateMasks,
                masksByTime,
                database,
            ),
        )
        t_spikes.start()

    # # ----------------------deadChannel-------------------------------------------
    if calcAll or calcAllDeadChan:
        # deadChannelMetric
        detrend = False  # defualt value for DC
        demean = True  # defualt value for DC
        taper = 0  # defualt value for DC
        st_dc = st.copy()
        t_dc = Process(
            target=_dead_channel_wrapper,
            args=(st_dc, network, station, database, logger, detrend, demean, taper),
        )
        t_dc.start()

        # deadChanADFMetric
        # t_dc.join()
        # st_dcADF = st.copy()
        # t_dc = Process(target=_dead_channelADF_wrapper,
        #                   args=(st_dcADF, network, station, logger, dcADF_win_size,
        #                         dcADF_pval_thresh, dcADF_threshold, generateMasks, masksByTime, dcADF_use_thresh,
        #                         database))
        # t_dc.start()

        # deadChanMeanMetric
        t_dc.join()
        st_dcM = st.copy()
        t_dc = Process(
            target=_dead_channelMean_wrapper,
            args=(
                st_dcM,
                network,
                station,
                logger,
                dcMean_win_size,
                dcMean_thresh,
                generateMasks,
                masksByTime,
                database,
            ),
        )
        t_dc.start()

    # # ----------------------calibration-------------------------------------------
    if calcAll or calcCal:
        st_cal = st.copy()
        _calibration_wrapper(
            st_cal, network, station, logger, cal_metric_store, database
        )

    # # ----------------------plots--------------------------------

    if plots:
        if calcAll:
            t_amp.join()
            t_basic.join()
            t_gap.join()
            t_snr.join()
            t_stalta.join()
            t_soh.join()
            t_dc.join()
            t_spikes.join()
            t_offset.join()

            try:
                if t_transfer.is_alive():
                    t_transfer.join()
            except UnboundLocalError:
                pass
            try:
                if t_cross.is_alive():
                    t_cross.join()
            except UnboundLocalError:
                pass
            try:
                if t_corr.is_alive():
                    t_corr.join()
            except UnboundLocalError:
                pass

        elif calcAll == None:
            pass
        else:
            if calcAmp:
                t_amp.join()
            if calcBasic:
                t_basic.join()
            if calcDcOffset:
                t_offset.join()
            if calcGap:
                t_gap.join()
            if calcSNR:
                t_snr.join()
            if calcStalta:
                t_stalta.join()
            if calcSOH:
                t_soh.join()
            if calcAllDeadChan:
                t_dc.join()
            if calcTransfer:
                t_transfer.join()
            if calcSpikes:
                t_spikes.join()

        print("-----------------------------------------------")
        print(("Plotting " + network + "." + station + " Jul Date: " + date))
        print("-----------------------------------------------")
        # st_plot = st.copy()

        dailyPdfplots(
            database,
            pdfModel,
            sta_dir + "/dailyPDFgrid",
            sta_dir + "/dailyPDFline",
            per_arr,
            network=network,
            station=station,
            channel=None,
            session=session,
        )
        print(("Finished PDFgrid and Line plots: " + network + "." + station))

        # st_noisePlot = st.copy()
        stationNoiseModel(
            database,
            plot=True,
            fname=sta_dir + "/stationNoiseModel.png",
            network=network,
            station=station,
            session=session,
        )
        print("Finished stationNoisePlot: " + network + "." + station)

        psdPlot(
            database,
            style="psd",
            f_name=sta_dir + "/psdPlot",
            showNoiseModel=showNoiseModel,
            showMaxMin=showMaxMin,
            showMode=showMode,
            showMean=showMean,
            showMedian=showMedian,
            showEnvelope=showEnvelope,
            envelopeType=envelopeType,
            showSingle=showSingle,
            singleType=singleType,
            ylo=ylo,
            yhi=yhi,
            pcolor=pqlx,
            timespan=timespan,
            network=network,
            station=station,
            channel=None,
            session=session,
        )
        print("Finished psdPlot: " + network + "." + station)

        psdPlot(
            database,
            "pdf",
            f_name=sta_dir + "/pdfPlot",
            showNoiseModel=False,
            showMaxMin=False,
            showMode=False,
            showMean=False,
            showMedian=False,
            showEnvelope=False,
            envelopeType=False,
            showSingle=False,
            singleType=False,
            ylo=ylo,
            yhi=yhi,
            pcolor=pqlx,
            timespan=timespan,
            network=network,
            station=station,
            channel=None,
            session=session,
        )
        print("Finished pdfPlot: " + network + "." + station)

        # Station ranking plots
        # try:
        #     if len(os.listdir(net_dir)) >= min_stations:
        #         if all([os.path.isdir(i) for i in os.listdir(net_dir)]):
        #             if all([len(i) == 3 or len(i) == 4 for i in os.listdir(net_dir)]):
        #                 pycheronStationRankingPlot(net_dir, pdfModel, rank_by)
        #                 print 'Finished Network Ranking Plot'
        # except ValueError:
        #     pass

        # moving files into correct folders
        # for fname in os.listdir(sta_dir):
        #     for j in un_chans:
        #         if fname.endswith(j + '.png'):
        #             shutil.move(sta_dir + '/' + fname, sta_dir + '/' + j + '/' + fname)

        print("-----------------------------------------------")
        print(("Finished Metric Calculations for " + network + "." + station))
        print("-----------------------------------------------")

    else:
        # if calcAll is false and calculating individual metrics, some processes may not be up
        if not calcAll:
            if calcAmp:
                t_amp.join()
            if calcBasic:
                t_basic.join()
            if calcDcOffset:
                t_offset.join()
            if calcGap:
                t_gap.join()
            if calcSNR:
                t_snr.join()
            if calcStalta:
                t_stalta.join()
            if calcSOH:
                t_soh.join()
            if calcAllDeadChan:
                t_dc.join()
            if calcTransfer:
                t_transfer.join()
            if calcSpikes:
                t_spikes.join()

        else:
            t_amp.join()
            t_basic.join()
            t_gap.join()
            t_snr.join()
            t_stalta.join()
            t_soh.join()
            t_dc.join()
            t_spikes.join()
            t_offset.join()
            t_transfer.join()
            t_cross.join()
            t_corr.join()

        # remove tmp MSEED file needed for SOH
        print("-----------------------------------------------")
        print(("Finished Metric Calculations for " + network + "." + station))
        print("-----------------------------------------------")


# -----------------------------------------------------------------------------
# Metric Wrapper Codes
# -----------------------------------------------------------------------------


def _psd_wrapper(
    st,
    expLoPeriod,
    expHiPeriod,
    linLoPeriod,
    linHiPeriod,
    evalresp,
    generateMasks,
    maskByTimes,
    dcExpThreshold,
    pctBelowNoiseThreshold,
    pctAboveNoiseThreshold,
    rmsThreshold,
    dcLinThreshold,
    num_gaps,
    pctBelowNoiseThresholdRESP,
    pctAboveNoiseThresholdRESP,
    processes,
    network,
    station,
    logger,
    byHourOn,
    database,
):
    """
    Wrapper for psdMetric

    If results are returned, 3 files (per channel) are created and saved to a Stream's channel dir. In addition, psds
    are saved in .npy files to be used later in plotting or other calculations, so psds wont have to be calculated again

    Output Metric filename: <network_dir>/<station_dir>/<channel_dir>/psdMetric_<channel>_<date>.csv
    Output Corrected PSD filename: <network_dir>/<station_dir>/<channel_dir>/corrected_psds_<channel>_<date>.csv
    Output PDF filename: <network_dir>/<station_dir>/<channel_dir>/pdfs_<channel>_<date>.csv
    Out PSD folder: <network_dir>/<station_dir>/PSDS
    """

    print("Calculating PSDS. This could take a while...")
    # calculating psdMetric
    psds = psdMetric(
        st,
        expLoPeriod=expLoPeriod,
        expHiPeriod=expHiPeriod,
        linLoPeriod=linLoPeriod,
        linHiPeriod=linHiPeriod,
        evalresp=evalresp,
        generateMasks=generateMasks,
        masksByTime=maskByTimes,
        dcExpThreshold=dcExpThreshold,
        pctBelowNoiseThreshold=pctBelowNoiseThreshold,
        pctAboveNoiseThreshold=pctAboveNoiseThreshold,
        rmsThreshold=rmsThreshold,
        dcLinThreshold=dcLinThreshold,
        num_gaps=num_gaps,
        pctBelowNoiseThresholdRESP=pctBelowNoiseThresholdRESP,
        pctAboveNoiseThresholdRESP=pctAboveNoiseThresholdRESP,
        processes=processes,
        logger=logger,
        byHourOn=byHourOn,
        database=database,
    )
    # if psds returned
    if psds:
        print(("callPycheron(): Finished psdMetric for: " + network + "." + station))
        logger.log("callPycheron(): Finished psdMetric for: " + network + "." + station)
    else:
        logger.warn("callPycheron(): No results returned for psdMetric()... Skipping")
        print("callPycheron(): No results returned for psdMetric()... Skipping")


def _basic_stats_wrapper(
    st,
    rmsThreshold,
    maxThreshold,
    minThreshold,
    medianThreshold,
    meanThreshold,
    varianceThreshold,
    stdThreshold,
    generateMasks,
    masksByTime,
    network,
    station,
    logger,
    database,
):
    """
    Wrapper code for basicStatsMetric

    If a result is returned from basicStatsMetric, a CSV file is written out to the Stream's channel folder

    Output filename: <network_dir>/<station_dir>/<channel_dir>/crossCorrMetric_<channel>_<date>.csv

    """
    # calculate basicStats
    basic_stats = basicStatsMetric(
        st,
        rmsThreshold=rmsThreshold,
        maxThreshold=maxThreshold,
        minThreshold=minThreshold,
        medianThreshold=medianThreshold,
        meanThreshold=meanThreshold,
        varianceThreshold=varianceThreshold,
        stdThreshold=stdThreshold,
        generateMasks=generateMasks,
        masksByTime=masksByTime,
        logger=logger,
        database=database,
    )

    if basic_stats:
        logger.log(
            "callPycheron(): Finished basicStatsMetric for: " + network + "." + station
        )
        print(
            (
                "callPycheron(): Finished basicStatsMetric for: "
                + network
                + "."
                + station
            )
        )
    else:
        # if basicStats does not return anything
        logger.warn(
            "callPycheron(): basicStatsMetric - No results returned... Skipping"
        )
        print("callPycheron(): basicStatsMetric - No results returned... Skipping")


def _corr_wrapper(stM, network, station, logger, database):
    """
    Wrapper for correlationMetric.

    If a result is returned from correlationMetric, a CSV file is written out to the Stream's station folder

    Output filename: <network_dir>/<station_dir>/correlationMetric_<date>.csv
    """
    # must have more than 2 traces in stream.
    if len(stM) >= 2:
        # calculating metric
        corrMetric = []
        # TODO: Make sure this matches all possible combinations, because stM might not be in order
        # loop through stream, get trace [i] and trace [i+1]
        for i in range(len(stM) - 1):
            tr1 = stM[i]
            tr2 = stM[i + 1]
            corrMetric = correlationMetric(tr1, tr2, logger, database)
            # if d returns result, add to list

        # if list contains results
        if corrMetric:
            logger.log(
                "callPycheron(): Finished correlationMetric for: "
                + network
                + "."
                + station
            )
            print(
                (
                    "callPycheron(): Finished correlationMetric for: "
                    + network
                    + "."
                    + station
                )
            )
        else:
            logger.warn(
                "callPycheron(): correlationMetric - No matching Traces... Skipping"
            )
            print("callPycheron(): correlationMetric - No matching Traces... Skipping")
    else:
        # if len is < 2
        logger.warn(
            "callPycheron(): correlationMetric - Not enought traces in Stream... Skipping"
        )
        print(
            "callPycheron(): correlationMetric - Not enought traces in Stream... Skipping"
        )


def _cross_corr_wrapper(
    stM,
    freqmin,
    maxLagSecs,
    filt,
    freqmax,
    corners,
    zerophase,
    maxorder,
    ba,
    freq_passband,
    network,
    station,
    logger,
    database,
):
    """
    Wrapper for crossCorrelationMetric.

    If a result is returned from crossCorrelationMetric, a CSV file is written out to the Stream's station folder

    Output filename: <network_dir>/<station_dir>/crossCorrMetric_<date>.csv
    """
    # must have two or more traces in stream
    if len(stM) >= 2:
        if freqmin is None:
            freqmin = 0.1
        # calculating metric
        crossCorr = []
        # TODO: Make sure this matches all possible combinations, because stM might not be in order
        # loop through stream, get trace [i] and trace [i+1]
        for i in range(len(stM) - 1):
            tr1 = stM[i]
            tr2 = stM[i + 1]
            # calc corssCorrMetric
            crossCorr = crossCorrMetric(
                tr1,
                tr2,
                maxLagSecs=maxLagSecs,
                filt=filt,
                freqmin=freqmin,
                freqmax=freqmax,
                corners=corners,
                zerophase=zerophase,
                maxorder=maxorder,
                ba=ba,
                freq_passband=freq_passband,
                logger=logger,
                database=database,
            )

        # if results in list
        if crossCorr:
            logger.log(
                "callPycheron(): Finished crossCorrMetric for: "
                + network
                + "."
                + station
            )
            print(
                (
                    "callPycheron(): Finished crossCorrMetric for: "
                    + network
                    + "."
                    + station
                )
            )
        else:
            # if results are empty, no matching traces to correlate
            logger.warn(
                "callPycheron(): crossCorrMetric - No matching traces... Skipping Metric"
            )
            print(
                "callPycheron(): crossCorrMetric - No matching traces... Skipping Metric"
            )
    else:
        # if len < 2, not enough traces for metric
        logger.warn("callPycheron(): crossCorrMetric Not enought traces in Stream")
        print("callPycheron(): crossCorrMetric Not enought traces in Stream")


def _offset_wrapper(
    st,
    generateMasks,
    masksByTime,
    windowSecs,
    incrementSecs,
    threshold,
    network,
    station,
    logger,
    database,
):
    """
    Wrapper for DCOffSetTimesMetric.

    If a result is returned from DCOffSetTimesMetric, a CSV file is written out to the Stream's respective channel
    folder.

    Output filename: <network_dir>/<station_dir>/<channel_dir>/DCOffSetTimesMetric_<channel>_<date>.csv
    """
    # calculating dc offsets
    offset_times = DCOffSetTimesMetric(
        st,
        generateMasks=generateMasks,
        masksByTime=masksByTime,
        windowSecs=windowSecs,
        incrementSecs=incrementSecs,
        threshold=threshold,
        logger=logger,
        database=database,
    )
    # if DCOffSetMetric returns results
    if offset_times:
        logger.log("callPycheron(): Finished DCOffset for:" + network + "." + station)
        print(("callPycheron(): Finished DCOffset for:" + network + "." + station))
    else:
        logger.warn(
            "callPycheron(): DCOffSetTimesMetric: No results returned... Skipping"
        )
        print("callPycheron(): DCOffSetTimesMetric: No results returned... Skipping")


def _gap_wrapper(
    st, generateMasks, separateMasks, masksByTime, network, station, logger, database
):
    """
    Wrapper for gapMetric

    If a result is returned from gapMetric, two CSV files are written. 1. A summary file for all the stations
    located in the Station folder, and 2. A values file, which contains the specific values and information about each
    gap, found in each respective channel folder. If there are not gaps a channel gap file will not be generated, only
    the summary file will be.

    Output summary filename: <station_directory>/gapsMetricSummary_<date>.csv
    Output values filename: <network_dir>/<station_dir>/<channel_dir>/gapsMetric_<channel>_<date>.csv
    """
    # calculate gaps metric.
    # gaps - list of individual gaps information
    # sum - summary information for station
    gapMetric(
        st,
        generateMasks=generateMasks,
        separateMasks=separateMasks,
        masksByTime=masksByTime,
        logger=logger,
        database=database,
    )

    logger.log("callPycheron(): Finished gapMetric for: " + network + "." + station)
    print(("callPycheron(): Finished gapMetric for: " + network + "." + station))


def _rep_amp_wrapper(
    st, minRep, generateMasks, network, station, logger, fortran, database
):
    """
    Wrapper for repeatedAmplitudeMetric

    If a result is returned from repeatedAmplitude Metric, two CSV files are written. 1. A summary file for all the stations
    located in the Station folder, and 2. A values file, which contains the specific values and information about each
    repeated amplitude, found in each respective channel folder. If there are no repeated amplitudes, a value file
    will not be generated

    Output summary filename: <station_directory>/repeatedAmplitudeSummary_<date>.csv
    Output values filename: <network_dir>/<station_dir>/<channel_dir>/repeatedAmplitudeValues_<channel>_<date>.csv
    """
    start = time.time()
    rep_amp = []
    # for each trace in stream, calculate repeatedAmplitudeMetric and append to list

    rep_amp = repeatedAmplitudeMetric(
        st,
        minRep=minRep,
        generateMasks=generateMasks,
        logger=logger,
        fortran=fortran,
        database=database,
    )

    # if results returned
    if rep_amp:
        logger.log(
            "callPycheron(): Finished repeatedAmplitude for: " + network + "." + station
        )
        print(
            (
                "callPycheron(): Finished repeatedAmplitude for: "
                + network
                + "."
                + station
            )
        )
    else:
        # if nothing returned
        logger.warn(
            "callPycheron(): repeatedAmplitudeMetric: No results returned... Skipping"
        )
        print(
            "callPycheron(): repeatedAmplitudeMetric: No results returned... Skipping"
        )
    timestart = (time.time() - start) / 60
    logger.log(
        "callPycheron(): repeatedAmplitudeMetric: Time in minutes: " + str(timestart)
    )


def _snr_wrapper(
    st,
    algorithmSNR,
    windowSecsSNR,
    snrThreshold,
    generateMasks,
    masksByTime,
    network,
    station,
    logger,
    database,
):
    """
    Wrapper for snrMetric.

    If a result is returned from snrMetric, a CSV file is written out to the Stream's respective channel
    folder.

    Output filename: <network_dir>/<station_dir>/<channel_dir>/snrMetric_<channel>_<date>.csv
    """
    # calculating snrMetric
    snr = snrMetric(
        st,
        algorithm=algorithmSNR,
        windowSecs=windowSecsSNR,
        snrThreshold=snrThreshold,
        generateMasks=generateMasks,
        masksByTime=masksByTime,
        logger=logger,
        database=database,
    )
    # if snrMetric returns a result, loop through unique channels and create/open a csv file to write out to.
    if snr:
        logger.log("callPycheron(): Finished snrMetric for: " + network + "." + station)
        print(("callPycheron(): Finished snrMetric for: " + network + "." + station))
    else:
        logger.warn("callPycheron(): snrMetric - No results returned... Skipping")
        # if nothing is returned
        print("callPycheron(): snrMetric - No results returned... Skipping")


def _soh_wrapper(
    st, data_quality, activity, io_clock, network, station, logger, database
):
    """
    Wrapper for sohMetric.

    If a result is returned from sohMetric, a CSV file is written out to the Stream's respective station
    folder.

    Output filename: <network_dir>/<station_dir>/sohMetric_<date>.csv
    """
    # calculate sohMetric
    soh = sohMetric(
        st,
        data_quality=data_quality,
        activity=activity,
        io_clock=io_clock,
        logger=logger,
        database=database,
    )
    # if soh returns result, create/open a csv file to write out to.
    if soh != None:
        logger.log("callPycheron(): Finished sohMetric for: " + network + "." + station)
        print(("callPycheron(): Finished sohMetric for: " + network + "." + station))

    else:
        logger.warn("callPycheron(): sohMetric returned None")
        print("callPycheron(): No results returned... Skipping")


def _stalta_wrapper(
    st,
    staSecs,
    ltaSecs,
    increment,
    algorithmSTA,
    network,
    station,
    logger,
    fortran,
    database,
):
    """
    Wrapper for staltaMetric.

    If a result is returned from staltaMetric, a CSV file is written out to the Stream's respective channel
    folder.

    Output filename: <network_dir>/<station_dir>/<channel_dir>/staltaMetric_<channel>_<date>.csv
    """
    start = time.time()
    # calculate stalta metric
    stalta = staltaMetric(
        st,
        staSecs=staSecs,
        ltaSecs=ltaSecs,
        increment=increment,
        algorithm=algorithmSTA,
        logger=logger,
        fortran=fortran,
        database=database,
    )
    # if staltaMetric returns a result, loop through unique channels and create/open a csv file to write out to.
    if stalta:
        logger.log(
            "callPycheron(): Finished staltaMetric for: " + network + "." + station
        )
        print(("callPycheron(): Finished staltaMetric for: " + network + "." + station))
    else:
        logger.warn("callPycheron(): staltaMetric: No results returned... Skipping")
        print("callPycheron(): staltaMetric: No results returned... Skipping")
    timestart = (time.time() - start) / 60
    logger.log("callPycheron(): staltaMetric: Time in minutes: " + str(timestart))


def _transfer_wrapper(st, network, station, logger, database):
    """
    Wrapper for transferFunctionMetric.

    If a stream have >= 2 traces, and transferFunction returns a results, a CSV file is written out to the Stream's
    respective channel folder

    Output filename: <network_dir>/<station_dir>/<channel_dir>/transferFunctionMetric_<channel>_<date>.csv
    """
    start = time.time()
    # Needs >= 2 traces in stream too compare for metric
    if len(st) >= 2:
        transMetric = []
        # Loop through stream, get stream [i] and stream  [i + 1]
        for i in range(len(st) - 1):
            tr1 = st[i]
            tr2 = st[i + 1]
            transMetric = transferFunctionMetric(tr1, tr2, logger, database=database)

        # If result is returned, loop through unique channels and create/open CSV file.
        if transMetric:
            logger.log(
                "callPycheron(): Finished transferFunctionMetric for: "
                + network
                + "."
                + station
            )
            print(
                (
                    "callPycheron(): Finished transferFunctionMetric for: "
                    + network
                    + "."
                    + station
                )
            )
        else:
            # if result returns nothing
            logger.warn(
                "callPycheron(): transferFunctionMetric: No Matching Traces... Skipping Metric"
            )
            print(
                "callPycheron(): transferFunctionMetric: No Matching Traces... Skipping Metric"
            )
    else:
        logger.warn(
            "callPycheron(): transferFunctionMetric: Not enought traces in Stream"
        )
        print("callPycheron(): transferFunctionMetric: Not enought traces in Stream")
    timestart = (time.time() - start) / 60
    logger.log(
        "callPycheron(): transferFunctionMetric: Time in minutes: " + str(timestart)
    )


def _spikes_wrapper(
    st,
    windowSize,
    thresholdSpikes,
    selectivity,
    fixedThreshold,
    processesSpikes,
    network,
    station,
    logger,
    fortran,
    generateMasks,
    masksByTime,
    database,
):
    """
    Wrapper for spikesMetric.

    If a result is returned from spikesMetric, a CSV file is written out to the Stream's respective channel
    folder.

    Output filename: <network_dir>/<station_dir>/<channel_dir>/spikesMetric_<channel>_<date>.csv
    """
    start = time.time()
    # calculate spikes metric
    spikes = spikesMetric(
        st,
        windowSize=windowSize,
        threshold=thresholdSpikes,
        selectivity=selectivity,
        fixedThreshold=fixedThreshold,
        processes=processesSpikes,
        logger=logger,
        fortran=fortran,
        generateMasks=generateMasks,
        masksByTime=masksByTime,
        database=database,
    )

    # if spikesMetric returns a result, loop through unique channels and create/open a csv file to write out to.
    if spikes:
        logger.log(
            "callPycheron(): Finished spikesMetric for: " + network + "." + station
        )
        print(("callPycheron(): Finished spikesMetric for: " + network + "." + station))
    else:
        # if nothing is returned from metric
        logger.warn("callPycheron(): spikesMetric: No result returned... Skipping")
        print("callPycheron(): spikesMetric: No result returned... Skipping")
    timestart = (time.time() - start) / 60
    logger.log("callPycheron(): spikesMetric: Time in minutes: " + str(timestart))


def _dead_channel_wrapper(
    st, network, station, database, logger=None, detrend=True, demean=True, taper=0.1
):
    """
    Wrapper for deadChannelMetric.

    If a result is returned from deadChannel Metric, a CSV file is written out to the Stream's respective channel
    folder.

    Output filename: deadChannelMetric_<channel>_<date>.csv
    """
    # calculate deadChannelMetric
    dc = deadChannelMetric(st, logger, database=database)
    # If result is returned, loop through unique channels and create/open CSV file.
    if dc:
        logger.log(
            "callPycheron(): Finished deadChannelMetric for: " + network + "." + station
        )
        print(
            (
                "callPycheron(): Finished deadChannelMetric for: "
                + network
                + "."
                + station
            )
        )
    else:
        # if no results are returned
        logger.warn(
            "callPycheron(): deadChannelMetric - No results returned... Skipping"
        )
        print("callPycheron(): deadChannelMetric - No results returned... Skipping")


def _dead_channelADF_wrapper(
    st,
    network,
    station,
    logger,
    win_size,
    pval_threshold,
    threshold,
    generateMasks,
    masksByTime,
    use_thesh,
    database,
):
    """
    Wrapper for deadChanADFMetric.

    If a result is returned from deadChanADFMetric, a CSV file is written out to the Stream's respective channel
    folder.

    Output filename: deadChanADFMetric_<channel>_<date>.csv

    """

    dc = deadChanADFMetric(
        st,
        win_size,
        pval_threshold,
        threshold,
        generateMasks,
        masksByTime,
        use_thesh,
        database=database,
    )

    if dc:
        logger.log(
            "callPycheron(): Finished deadChanADFMetric for: " + network + "." + station
        )
        print(
            (
                "callPycheron(): Finished deadChanADFMetric for: "
                + network
                + "."
                + station
            )
        )
    else:
        # if no results are returned
        logger.warn(
            "callPycheron(): deadChanADFMetric - No results returned... Skipping"
        )
        print("callPycheron(): deadChanADFMetric - No results returned... Skipping")


def _dead_channelMean_wrapper(
    st,
    network,
    station,
    logger,
    win_size,
    threshold,
    generateMasks,
    masksByTime,
    database,
):
    """
    Wrapper for deadChanADFMetric.

    If a result is returned from deadChanADFMetric, a CSV file is written out to the Stream's respective channel
    folder.

    Output filename: deadChanMeanMetric_<channel>_<date>.csv

    """

    dc = deadChanMean(
        st, win_size, threshold, generateMasks, masksByTime, database=database
    )

    if dc:
        logger.log(
            "callPycheron(): Finished deadChanMeanMetric for: "
            + network
            + "."
            + station
        )
        print(
            (
                "callPycheron(): Finished deadChanMeanMetric for: "
                + network
                + "."
                + station
            )
        )
    else:
        # if no results are returned
        logger.warn(
            "callPycheron(): deadChanMeanMetric - No results returned... Skipping"
        )
        print("callPycheron(): deadChanMeanMetric - No results returned... Skipping")


def _calibration_wrapper(
    st, network, station, logger, metric_store=None, database=None
):
    """
    Wrapper for calibration metric

    If a result is returned from calibration, a CSV file is written out to the Stream's respective channel
    folder.

    Output filename: calibrationMetric_<channel>_<date>.csv
    """
    cal = calibrationMetric(st, metric_store=metric_store, database=database)

    if cal:
        logger.log(
            "callPycheron(): Finished calibrationMetric for: " + network + "." + station
        )
        print(
            (
                "callPycheron(): Finished calibrationMetric for: "
                + network
                + "."
                + station
            )
        )
    else:
        # if no results are returned
        logger.warn(
            "callPycheron(): Finished calibrationMetric for: " + network + "." + station
        )
        print(
            (
                "callPycheron(): Finished calibrationMetric for: "
                + network
                + "."
                + station
            )
        )


# -----------------------------------------------------------------------------
# Main Code
# -----------------------------------------------------------------------------


def main():
    import time
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description="callPycheronMetric")
    parser.add_argument(
        "config_file",
        type=str,
        help="Configuration file containing callPycheronMetric params",
    )

    args = parser.parse_args()

    if args:
        config = args.config_file

        if config.split(".")[-1] == "yaml":
            with open(config, "r") as ymlfile:
                cfg = yaml.load(ymlfile)
        else:
            print("Error: config file must be yaml")
            return
    else:
        print("Error: No config file")
        return

    for key, val in cfg.items():
        if val == "None" and key != "session":
            cfg[key] = None

    start = time.time()

    callPycheron(
        output_dir=cfg["output_dir"],
        data=cfg["data"],
        datatype=cfg["datatype"],
        calcAll=cfg["calcAll"],
        calcPsds=cfg["calcPsds"],
        calcBasic=cfg["calcBasic"],
        calcCorr=cfg["calcCorr"],
        calcCrossCorr=cfg["calcCrossCorr"],
        calcGap=cfg["calcGap"],
        calcAmp=cfg["calcAmp"],
        calcSNR=cfg["calcSNR"],
        calcSOH=cfg["calcSOH"],
        calcStalta=cfg["calcStalta"],
        calcDcOffset=cfg["calcDcOffset"],
        calcSpikes=cfg["calcSpikes"],
        calcAllDeadChan=cfg["calcAllDeadChan"],
        calcTransfer=cfg["calcTransfer"],
        calcCal=cfg["calcCal"],
        network=cfg["network"],
        station=cfg["station"],
        byDay=cfg["byDay"],
        startdate=cfg["startdate"],
        enddate=cfg["enddate"],
        jul_start=cfg["jul_start"],
        jul_end=cfg["jul_end"],
        generateMasks=cfg["generateMasks"],
        masksByTime=cfg["masksByTime"],
        rmsThreshold=cfg["rmsThreshold"],
        maxThreshold=cfg["maxThreshold"],
        minThreshold=cfg["minThreshold"],
        medianThreshold=cfg["medianThreshold"],
        meanThreshold=cfg["meanThreshold"],
        varianceThreshold=cfg["varianceThreshold"],
        stdThreshold=cfg["stdThreshold"],
        maxLagSecs=cfg["maxLagSecs"],
        filt=cfg["filt"],
        freqmin=cfg["freqmin"],
        freqmax=cfg["freqmax"],
        corners=cfg["corners"],
        zerophase=cfg["zerophase"],
        maxorder=cfg["maxorder"],
        ba=cfg["ba"],
        freq_passband=cfg["freq_passband"],
        windowSecs=cfg["windowSecs"],
        incrementSecs=cfg["incrementSecs"],
        threshold=cfg["threshold"],
        separateMasks=cfg["separateMasks"],
        completeDay=cfg["completeDay"],
        expLoPeriod=cfg["expLoPeriod"],
        expHiPeriod=cfg["expHiPeriod"],
        linLoPeriod=cfg["linLoPeriod"],
        linHiPeriod=cfg["linHiPeriod"],
        evalresp=cfg["evalresp"],
        dcExpThreshold=cfg["dcExpThreshold"],
        pctBelowNoiseThreshold=cfg["pctBelowNoiseThreshold"],
        pctAboveNoiseThreshold=cfg["pctAboveNoiseThreshold"],
        dcLinThreshold=cfg["dcLinThreshold"],
        num_gaps=cfg["num_gaps"],
        pctBelowNoiseThresholdRESP=cfg["pctBelowNoiseThresholdRESP"],
        pctAboveNoiseThresholdRESP=cfg["pctAboveNoiseThresholdRESP"],
        minRep=cfg["minRep"],
        algorithmSNR=cfg["algorithmSNR"],
        windowSecsSNR=cfg["windowSecsSNR"],
        snrThreshold=cfg["snrThreshold"],
        data_quality=cfg["data_quality"],
        activity=cfg["activity"],
        io_clock=cfg["io_clock"],
        windowSize=cfg["windowSize"],
        thresholdSpikes=cfg["thresholdSpikes"],
        selectivity=cfg["selectivity"],
        fixedThreshold=cfg["fixedThreshold"],
        staSecs=cfg["staSecs"],
        ltaSecs=cfg["ltaSecs"],
        increment=cfg["increment"],
        algorithmSTA=cfg["algorithmSTA"],
        plots=cfg["plots"],
        pdfModel=cfg["pdfModel"],
        per_arr=cfg["per_arr"],
        showNoiseModel=cfg["showNoiseModel"],
        showMaxMin=cfg["showMaxMin"],
        showMode=cfg["showMode"],
        showMean=cfg["showMean"],
        showMedian=cfg["showMedian"],
        showEnvelope=cfg["showEnvelope"],
        envelopeType=cfg["envelopeType"],
        showSingle=cfg["showSingle"],
        singleType=cfg["singleType"],
        ylo=cfg["ylo"],
        yhi=cfg["yhi"],
        min_stations=cfg["min_stations"],
        rank_by=cfg["rank_by"],
        processesPSD=cfg["processesPSD"],
        processesSpikes=cfg["processesSpikes"],
        log=cfg["log"],
        fortran=cfg["fortran"],
        timespan=cfg["timespan"],
        dcADF_win_size=cfg["dcADF_win_size"],
        dcADF_pval_thresh=cfg["dcADF_pval_thresh"],
        dcADF_threshold=cfg["dcADF_threshold"],
        dcADF_use_thresh=cfg["dcADF_use_thresh"],
        dcMean_win_size=cfg["dcMean_win_size"],
        dcMean_thresh=cfg["dcMean_thresh"],
        cal_metric_store=cfg["cal_metric_store"],
        dcExpThresholdHour=cfg["dcExpThresholdHour"],
        dcLinThresholdHour=cfg["dcLinThresholdHour"],
        byHourOn=cfg["byHourOn"],
        database=cfg["database"],
        session=cfg["session"],
        overwrite=cfg["overwrite"],
        to_csv=cfg["to_csv"],
        stationStartAt=cfg["stationStartAt"],
    )

    timestart = (time.time() - start) / 60
    print(("Time in minutes: " + str(timestart)))


if __name__ == "__main__":
    main()
