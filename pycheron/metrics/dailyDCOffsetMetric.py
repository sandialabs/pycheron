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
from scipy.stats.mstats import mquantiles


def dailyDCOffSetMetric(
    d, offsetDays=5, outlierWindow=7, outlierThreshold=6.0, OutputType=1
):
    """
    Metric to process a dataframe with daily means from basicStatsMetric and return a vector of daily likelihoods
    that a DC shift occurred, e.g., identifies days with a jump in the signal mean.

    :param d: Dictionary of lists (could be fed back as a dataframe, depends on firebird team)
    :type d: dict
    :param offsetDays: Number of days used in calculating weighting functions
    :type offsetDays: int
    :param outlierWindow: Window size passed to findOutliers() function
    :type outlierWindow: int
    :param outlierThreshold:  Detection threshold passed to findOutliers() function
    :type outlierThreshold: float
    :param OutputType: if 1, return last day of valid values (`index=len(index-np.floor(outlierWindow/2))`); if 0, return all valid values (`indices=max(offsetDays,np.floor(outlierWindow/2))`):
                       `len(index)-np.floor(outlierWindow/2)`
    :type OutputType: bool

    :return: A list or list of dictionaries is returned for the last day `np.floor(outlierWindow/2)` (default 3rd day from last) in incoming data if outputType = 1
             (one list element), otherwise the first + offsetDays to last day `np.floor(outlierWindow/2)` multiple list elements, (one per day) is returned if outputType = 0.
    :rtype: list

    **Examples**

    .. code-block:: python

        # Example is just an example input, as this generated input doesn't exist, would come from a dataframe built up over time
        df = pd.read_csv('docs/test_data/test.csv')
        offsetDays = 5
        outlierWindow = 7
        outlierThreshold = 6
        OutputType = 0
        # Not shown, output type = 0 provides all days of incoming data (see out information above)
        dailyDCOffSetMetric (df, offsetDays,outlierWindow, outlierThreshold,OutputType)
        # output type = 1 provides last day of incoming data (see out information above) not shown
        # Not shown, output type = 0 provides all days of incoming data (see out information above). Would simply be last entry in outputType = 0.


    .. note::

        Prefer 60+ days of mean values to get a good estimate of the long term mean std. dev. After initial testing on
        stations in the IU network, a metric value > 10 appears to be indicative of a DC offset shift
        (this may vary across stations or networks and larger values may be preferred as indications of a
        potential station issue).

    """

    # # Convert dictionary of lists from basicStatsMetric to a dataframe
    df = pd.DataFrame(d)

    # Test whether dataFrame empty and if does not contain mean variable
    if df.empty:
        return
    if "mean" not in df:
        return

    a = np.asarray(df["mean"].copy())
    outliers = findOutliers(a, outlierWindow, outlierThreshold)
    cleanMean = np.asarray(a.copy())
    if outliers is not None:
        for i in outliers:
            cleanMean[i] = roll_median(a, outlierWindow, 1)[i]
    cleanMean = cleanMean[0 : int((len(cleanMean) - np.floor(outlierWindow / 2)))]

    metric = np.ones(len(cleanMean))

    # Have a minimum value to prevent occasional zeros associated with different lags from completely wiping out
    # large values
    diffMin = np.repeat(0.001, len(cleanMean))

    # Create vectors of daily differences with N=offsetDays increasing lags, multiplying them together and then
    # taking the N'th root
    dailyDiff = np.array([])

    # Create a daily metric from the lagged data with NaNs at the beginning.  Each date has the difference between that
    # date and the value 'i' days earlier
    for i in range(1, offsetDays + 1):
        dailyDiff = np.append(np.full(i, np.nan), abs(cleanMean[i:] - cleanMean[:-i]))

        try:
            # Multiplying them together weights those shifts over the previous N days
            metric *= np.maximum(diffMin, dailyDiff)
        except ValueError:
            print("dailyDCOffSetMetric: Not enough data to process")

    metric **= 1 / offsetDays

    # Scale the metric by the median of the rolling sd with a window size of offsetDays
    scaling = roll_sd(cleanMean, offsetDays, 1, "center")
    scal = mquantiles(scaling[~np.isnan(scaling)], 0.5)
    metric /= scal

    # Create list of metric values
    metricList = []

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
                    "value": metric[i],
                }
            metricList = metricL
        else:
            if np.isnan(metric[i]):
                continue
            else:
                metricL = {
                    "snclq": df["snclq"][i],
                    "start_time": df["start_time"][i],
                    "end_time": df["end_time"][i],
                    "metric_name": "dailyDCOffset",
                    "value": metric[i],
                }
            metricList.append(metricL)

    return metricList
