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

__all__ = ["dailyDCOffSetMetric"]

import numpy as np
import pandas as pd
from pycheron.psd.noise.findOutliers import findOutliers
from pycheron.rollseis.roll_median import roll_median
from pycheron.rollseis.roll_sd import roll_sd
from pycheron.util.logger import Logger
from scipy.stats.mstats import mquantiles
from pycheron.db.sqllite_db import Database


def dailyDCOffSetMetric(
    d,
    offsetDays=5,
    outlierWindow=7,
    outlierThreshold=6.0,
    OutputType=1,
    logger=None,
    database_config=None,
):
    """
    Metric to process a list of dictionaries with daily means from basicStatsMetric and return a vector of daily
    likelihoods that a DC shift occurred, i.e., identifies days with a jump in the signal mean.

    :param d: List of dictionaries
    :type d: list
    :param offsetDays: Number of days used in calculating weighting functions
    :type offsetDays: int
    :param outlierWindow: Window size passed to pycheron.psd.noise.findOutliers function
    :type outlierWindow: int
    :param outlierThreshold: Detection threshold passed to pycheron.psd.noise.findOutliers function
    :type outlierThreshold: float
    :param OutputType: if 1, return last day of valid values (`index=len(index-np.floor(outlierWindow/2))`); if 0,
                       return all valid values (`indices=max(offsetDays,np.floor(outlierWindow/2))`):
                       `len(index)-np.floor(outlierWindow/2)`
    :type OutputType: bool
    :param logger: If using a logger, (you must create one using the util.logger class)
    :type logger: pycheron.util.logger.Logger
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database

    :return: A list or list of dictionaries is returned for the `last day - np.floor(outlierWindow/2)`
             (default 3rd day from last) in incoming data if outputType = 1
             (one list element), otherwise the first + offsetDays to last day `np.floor(outlierWindow/2)` multiple list
             elements, (one per day) is returned if outputType = 0.
    :rtype: list

    * Code originally ported from IRISMustangMetrics R Cran Package
    (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, March 20). CRAN-Package IRISMustangMetrics.
    The Comprehensive R Archive Network. Retrieved from
    https://cran.r-project.org/web/packages/IRISMustangMetrics/index.html) and augmented and adapted for use within
    Pycheron

     ** Algorithm steps:**
      #. Convert input to a dataframe
      #. Perform basic integrity checks
      #. Find any outliers present and replace them with a rolling median value
      #. Calculate lagged differences in daily means
      #. Multiply lags together (metric *= np.max(diffMin, dailydiff) to weight shifts over the previous N days,
        then take Nth root (metric ** 1/offset days)
      #. Scale the metric by the median of the rolling standard deviation with a window size of offsetDays

      Calculates lagged differences in daily mean taken from input basicStatsMetric over the offsetDays supplied
    **Examples**

    .. code-block:: python

        # Example is just an example input, as this generated input doesn't exist, would come from a dataframe built up
        over time
        df = pd.read_csv('docs/test_data/test.csv')
        offsetDays = 5
        outlierWindow = 7
        outlierThreshold = 6
        OutputType = 0
        # Not shown, output type = 0 provides all days of incoming data (see out information above)
        dailyDCOffSetMetric (df, offsetDays, outlierWindow, outlierThreshold,OutputType)
        # output type = 1 provides last day of incoming data (see out information above) not shown
        # Not shown, output type = 0 provides all days of incoming data (see out information above). Would simply be
          last entry in outputType = 0.


    .. note::

        Notes from IRISMustangMetrics package
        (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf):
        "Prefer 60+ days of mean values to get a good estimate of the long term mean std. dev. After initial testing on
        stations in the IU network, a metric value > 10 appears to be indicative of a DC offset shift
        (this may vary across stations or networks and larger values may be preferred as indications of a
        potential station issue)."

    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # # Convert dictionary of lists from basicStatsMetric to a dataframe
    df = pd.DataFrame(d)

    # Test whether dataFrame is empty
    if df.empty:
        logger.error("dailyDCOffsetMetric(): dataframe is empty")
        return

    # Check if mean column exists in dataframe
    if "mean" not in df:
        logger.error("dailyDCOffsetMetric(): mean column does not exist in dataframe")
        return

    # Check if start_time column exists in dataframe
    if "start_time" not in df:
        logger.error("dailyDCOffsetMetric(): start_time column does not exist in dataframe")
        return

    # Check if end_time column exists in dataframe
    if "end_time" not in df:
        logger.error("dailyDCOffsetMetric(): end_time column does not exist in dataframe")
        return

    # Make an array copy of the mean column
    a = np.asarray(df["mean"].copy())

    # Calculate outliers
    outliers = findOutliers(a, outlierWindow, outlierThreshold)
    # Create a copy of the a array, which is a copy of the df's mean column
    cleanMean = np.asarray(a.copy())
    # If outliers are found, replace outlier values with rolling median values
    if outliers is not None:
        for i in outliers:
            cleanMean[i] = roll_median(a, outlierWindow, 1)[i]
    # The last np.floor(outlierWindow/2) are not checked in findOutliers, so remove
    cleanMean = cleanMean[0 : int((len(cleanMean) - np.floor(outlierWindow / 2)))]

    # Create an array full of 1's based on the length of cleanMean
    metric = np.ones(len(cleanMean))

    # Create minimum value array that has min value repeated by length of cleanMean. This is to prevent occasional zeros
    # associated with different lags from completely wiping out large values
    diffMin = np.repeat(0.001, len(cleanMean))

    # Create a daily diff from the lagged data with NaNs at the beginning.  Each date has the difference between that
    # date and the value 'i' days earlier.
    # Daily diff array is appended to as cycle through offsetDays
    for i in range(1, offsetDays + 1):
        dailyDiff = np.append(np.full(i, np.nan), abs(cleanMean[i:] - cleanMean[:-i]))

        try:
            # Multiply metric by np.max(diffMin, dailydiff) to weight shifts over the previous N days
            metric *= np.maximum(diffMin, dailyDiff)
        except ValueError:
            print("dailyDCOffSetMetric: Not enough data to process")

    # Take 1/offsetDays root of metric (i.e., Nth root)
    metric **= 1 / offsetDays

    # Scale metric by the median of the rolling sd. Window size is offsetDays, starting at index 1 and using
    # center alignment
    scaling = roll_sd(cleanMean, offsetDays, 1, "center")
    scal = mquantiles(scaling[~np.isnan(scaling)], 0.5)
    metric /= scal

    # Create list of metric values by looping through metric
    metricList = []

    # If OutputType = 1, then return last valid day of metric
    for i in range(len(metric)):
        if OutputType == 1:
            if np.isnan(metric[i]):
                continue
            else:
                metricL = {
                    "snclq": df["snclq"][i],
                    "start_time": df["start_time"][i],
                    "end_time": df["end_time"][i],
                    "metric_name": "dailyDCOffset",
                    "daily_dc_offset_value": metric[i],
                }
            metricList = metricL
        # If OutputType = 0, return all metrics
        else:
            if np.isnan(metric[i]):
                continue
            else:
                metricL = {
                    "snclq": df["snclq"][i],
                    "start_time": df["start_time"][i],
                    "end_time": df["end_time"][i],
                    "metric_name": "dailyDCOffset",
                    "daily_dc_offset_value": metric[i],
                }
            metricList.append(metricL)

    if not metricList:
        snclqs = [x for x in df["snclq"]]
        for i in range(len(snclqs)):
            metricList.append(
                {
                    "snclq": df["snclq"][i],
                    "start_time": df["start_time"][i],
                    "end_time": df["end_time"][i],
                    "metric_name": "dailyDCOffset",
                    "daily_dc_offset_value": None,
                }
            )

    # If database defined, insert metric information
    if database_config is not None:
        database = Database(**database_config)
        database.insert_metric(metricList)

    return metricList
