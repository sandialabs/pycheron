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

import numpy as np
from statsmodels.tsa.stattools import adfuller
from pycheron.util.masks import consecutive
from pycheron.db.sqllite_db import Database

__all__ = ["deadChanADFMetric"]


def deadChanADFMetric(
    st,
    win_size=500,
    pvalue_thresh=0.01,
    threshold=1.5,
    generateMasks=False,
    masksByTime=True,
    use_thresh="pvalue",
    database_config=None,
):
    """
    DeadChannelMetric test using the Augmented Dickey Fuller Test

    :param st: obspy stream object to run Augmented Dickey Fuller (ADF) Test on
    :type st: obspy.core.stream.Stream
    :param win_size: window size to increment the data by in samples, each window will have a ADF run on it
                     (assumes 100sps). If None, window size determined by sample rate of each trace in the stream object
    :type win_size: int
    :param pvalue_threshold: p value threshold to trigger ADF rejection of null hypothesis (e.g., time series is assumed
                             stationary with rejection of null hypothesis)
    :type pvalue_thresh: float
    :param threshold: threshold to trigger ADF rejection of null hypothesis using difference between DF test statistic
                      and 1% critical value
    :type threshold: float
    :param generateMasks: Generate boolean qc mask array
    :type generateMasks: bool
    :param masksByTime: Generate a time-based qc mask list
    :type masksByTime: bool
    :param use_thresh: String option for choosing whether to use the pvalue_threshold or the threshold for determining
                       whether null hypothesis fails. Options are `pvalue` or `diffs`.
    :type use_thresh: str
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database


    :return: Returns list of dictionaries that includes the following for each trace within the stream:

            * dead_channel_times (`list`) - Lists of dicts with dead channel start/endtime
            * end_time (`str`) - end time of trace
            * start_time (`str`) - start time of trace
            * masks (`list` or `numpy.ndarray`) - masks for that trace, which is either a boolean array or a list of
              start/stop times in UTC
            * metric_name (`str`) - deadChanADFMetric string
            * snlcq (`str`) - station network location channel information for the specific trace

    :rtype: list

    The Augmented Dickey Fuller (ADF) Test is unit root test for stationarity. Rejection of the ADF Test implies that
    the data within the window are stationary. Acceptance of the ADF test, implies unit root, or a non-stationary
    time series. The ADF in this case is used to find areas within the time series that are stationary (statistical
    properties such as mean, variance, autocorrelation, etc. are all constant over time) and mark them as a deadChannel
    QC issue.

    ** Resources for ADF Method **:

    #. Dickey, D.A., and W.A. Fuller (1979). Distribution of the estimators for autoregressive time series
       with a unit root, J Am Stat Assoc. 74(366a) 427-431.

    #. MacKinnon, J.G. (1994). Approximate asymptotic distribution functions for unit-root and cointegration tests,
       J Bus Econ Stat. 12(2) 167-176.

    #. MacKinnon, J.G. (2010). Critical values for cointegration tests (No. 1227), Queen's Economics Department Working
       Paper.




    **Examples**

    .. code-block:: python


    """

    # Loop through stream object and calculate a moving average for each trace based on user defined window size and
    # period -- Using pandas built in moving average because very fast

    # initialize out list
    out = []
    for s in range(len(st)):

        # Pull out trace one by one and loop through
        tr = st[s]

        # Obtain appropriate window size based on the sample rate
        if win_size is None:
            if int(tr.stats.sampling_rate) == 100:
                win_size = 500
            elif int(tr.stats.sampling_rate) == 40:
                win_size = 400
            elif int(tr.stats.sampling_rate) == 20:
                win_size = 200
            # TODO: May want to think about changing this size for other sample rates -- maybe split into higher than
            # TODO: 100sps and lower than 20sps?
            else:
                win_size = 500

        # Initialize empty list to append means to
        adfuller_rolling = []

        # Loop through data and calculate the mean of each window
        for n in np.arange(0, len(tr.data) - win_size, win_size):
            # Define n_start, n_end
            n_start = n
            n_end = n + win_size

            # Calculate rolling augmented dickey fuller
            adfuller_rolling.append(adfuller(tr.data[n_start:n_end]))

        # Convert list to array
        adfuller_rolling = np.asarray(adfuller_rolling)

        # Get differences between adf test statistics and the 1% critical value for every element
        diffs = np.zeros(len(adfuller_rolling))
        p_inds = np.zeros(len(adfuller_rolling))
        for adf in range(len(adfuller_rolling)):
            diffs[adf] = adfuller_rolling[adf][0] - adfuller_rolling[adf][4]["1%"]
            p_inds[adf] = abs(adfuller_rolling[adf][1]) <= pvalue_thresh

        # Find indices where differences are greater than threshold
        inds = np.argwhere(diffs <= -threshold).ravel()

        # If no inds are returned, skip
        if inds.size != 0:

            # Check p-values for arrays
            p_inds = np.argwhere(p_inds == 1).ravel()

            # Initialize deadMasks to None in the case where user turns off masking
            deadMask = None

            # Convert indices to time for masks.
            # If consecutive indices, use start/end time. If only one use ind as start
            # and end time. E.g., for long periods of time use first index and last index until a break in being below
            # the threshold
            # Initialize dead channel times list
            dead_times = []

            # Group array into consecutive element arrays to identify stop/start times easier
            if use_thresh == "pvalue":
                cons = consecutive(p_inds)
            else:
                cons = consecutive(inds)

            # Loop through range of consecutive elements to get the start/end time of dead Channel issues
            for c in range(len(cons)):
                # Skip filter flukes and only process consecutive elements greater than con_thresh due to Type1
                # errors in ADF test
                if int(tr.stats.sampling_rate) == 100:
                    # window filtering threshold is ~2.1 min to eliminate flukes
                    # TODO: 10 min is 120 for win_size 500, 100sps might be more effective to window out issues
                    con_thresh = (2.1 * 60 * 100) / win_size
                elif int(tr.stats.sampling_rate) == 40:
                    con_thresh = (2.1 * 60 * 40) / win_size
                elif int(tr.stats.sampling_rate) == 20:
                    con_thresh = (2.1 * 60 * 20) / win_size
                else:
                    con_thresh = (2.1 * 60 * tr.stats.sampling_rate) / win_size
                if len(cons[c]) >= con_thresh:
                    start = tr.stats.starttime + (cons[c][0] * win_size) * 1.0 / tr.stats.sampling_rate
                    end = tr.stats.starttime + (cons[c][-1] * win_size) * 1.0 / tr.stats.sampling_rate

                    # Fill dictionary for each of these start/time intervals
                    d = {
                        "dead_starttime": start.isoformat(),
                        "dead_endtime": end.isoformat(),
                    }
                    dead_times.append(d)

            # If generateMasks set to True, boolean masks will be created for QC issues -- skip over single elements of
            # array as these are likely flukes
            if generateMasks:
                # If masksByTime set to True, masks with start/stop times will be created for QC issues
                if masksByTime:
                    # Create empty list to append to
                    deadMask = []
                    # Loop through consecutive element blocks, skip single element issues.
                    # Append Starttime and Endtime to list, then convert to array
                    for c in range(len(cons)):
                        if int(tr.stats.sampling_rate) == 100:
                            # window filtering threshold is ~2.1 min to eliminate flukes
                            # potential TODO: 10 min is 120 for win_size 500,
                            # 100sps might be more effective to window out issues
                            con_thresh = (2.1 * 60 * 100) / win_size
                        elif int(tr.stats.sampling_rate) == 40:
                            con_thresh = (2.1 * 60 * 40) / win_size
                        elif int(tr.stats.sampling_rate) == 20:
                            con_thresh = (2.1 * 60 * 20) / win_size
                        else:
                            con_thresh = (2.1 * 60 * tr.stats.sampling_rate) / win_size
                        if len(cons[c]) >= con_thresh:
                            s = tr.stats.starttime + (cons[c][0] * win_size) * 1.0 / tr.stats.sampling_rate
                            e = tr.stats.starttime + cons[c][-1] * win_size * 1.0 / tr.stats.sampling_rate
                            deadMask.append({"Startime": s.isoformat(), "Endtime": e.isoformat()})

                        # Append the general information below
                        dt = {
                            "snclq": tr.get_id(),
                            "start_time": tr.stats.starttime.isoformat(),
                            "end_time": tr.stats.endtime.isoformat(),
                            "dead_channel_times": dead_times if dead_times else None,
                            "masks": deadMask,
                            "metric_name": "deadChanADFMetric",
                        }
                        out.append(dt)
                else:
                    # if maskByTime set to false, create binary array
                    deadMask = np.zeros(len(tr.data), dtype=int)
                    for c in range(len(cons)):
                        if int(tr.stats.sampling_rate) == 100:
                            # window filtering threshold is ~2.1 min to eliminate flukes
                            # potential TODO: 10 min is 120 for win_size 500,
                            # 100sps might be more effective to window out issues
                            con_thresh = (2.1 * 60 * 100) / win_size
                        elif int(tr.stats.sampling_rate) == 40:
                            con_thresh = (2.1 * 60 * 40) / win_size
                        elif int(tr.stats.sampling_rate) == 20:
                            con_thresh = (2.1 * 60 * 20) / win_size
                        else:
                            con_thresh = (2.1 * 60 * tr.stats.sampling_rate) / win_size
                        if len(cons[c]) >= con_thresh:
                            deadMask[cons[c] * win_size] = 1

                    # Append the general information below
                    dt = {
                        "snclq": tr.get_id(),
                        "start_time": tr.stats.starttime.isoformat(),
                        "end_time": tr.stats.endtime.isoformat(),
                        "dead_channel_times": dead_times if dead_times else None,
                        "masks": deadMask,
                        "metric_name": "deadChanADFMetric",
                    }
                    out.append(dt)

            else:

                # Append the general information below for every trace, as it will be the same
                dt = {
                    "snclq": tr.get_id(),
                    "start_time": tr.stats.starttime.isoformat(),
                    "end_time": tr.stats.endtime.isoformat(),
                    "dead_channel_times": dead_times if dead_times else None,
                    "masks": "No Masks Created",
                    "metric_name": "deadChanADFMetric",
                }
                out.append(dt)
        else:
            dt = {
                "snclq": tr.get_id(),
                "start_time": tr.stats.starttime.isoformat(),
                "end_time": tr.stats.endtime.isoformat(),
                "dead_channel_times": None,
                "masks": "No Masks Created",
                "metric_name": "deadChanADFMetric",
            }
            out.append(dt)

    # If database defined, insert metric information
    if database_config is not None:
        database = Database(**database_config)
        database.insert_metric(out)

    return out
