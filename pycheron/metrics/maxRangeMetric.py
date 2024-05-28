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

__all__ = ["maxRange"]

import numpy as np
from pycheron.util.logger import Logger
from pycheron.rollseis.roll_range import roll_range
from pycheron.db.sqllite_db import Database


def maxRange(st, window=300, increment=150, database_config=None, logger=None):
    """
    Calculates the maximum sample range (difference between the largest and smallest sample value)
    for the given stream object over a user defined rolling window, incremented by a user defined increment

    :param st: ObsPy Stream object
    :type st: obspy.core.stream.Stream
    :param window: length of rolling window (in s) to evaluate min/max sample values (default = 300s)
    :type window: integer
    :param increment: number of seconds to advance window for each calculation (default = 150s)
    :type increment: integer
    :param database_config: dictionary containing the necessary parameters to create
                            a pycheron Database object. 
                            These include "db_name", "session_name", "overwrite", "manual", "wfdb_conn"
    :type database_config: dict
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger

    :return: list of dictionaries of the largest max range encountered for each trace object within the stream object
             Calculates the difference between the largest and smallest samples within a rolling window incrementing
             through the time series.
             The following keys exist:

            * snclq (`str`): station network channel location quality code information for specfic trace
            * start_time (`str`): start time of trace
            * end_time (`str`): end time of trace
            * maxRange (`float`): largest max range encountered for the respective trace object
            * metric_name (`str`): maxRangeMetric string

    :rtype: list

    * Code originally ported from IRISMustangMetrics R Cran Package
    (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, March 20). CRAN-Package IRISMustangMetrics.
    The Comprehensive R Archive Network. Retrieved from
    https://cran.r-project.org/web/packages/IRISMustangMetrics/index.html) and augmented and adapted for use within
    Pycheron
    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Set up empty list for appending metric information for each trace
    d = []

    # Warn user if they are selecting sub-optimal increment value, outside of given best practices
    if increment >= window:
        logger.warn("For optimal results, increment value should be less than window value")

    # Merge any traces within the stream object that may share the same snclq
    st.merge()

    # Loop through all traces in provided stream object
    for tr in st:
        # Get snclq, starttime, endtime information for each trace
        snclq = tr.get_id()
        starttime = tr.stats.starttime
        endtime = tr.stats.endtime

        # Get integer number of samples and calculate increment needed based on sampling rate
        n_samps = np.round(window * tr.stats.sampling_rate)
        n_incr = np.round(increment * tr.stats.sampling_rate)

        # Ensure that when creating the last rolling segment that it breaks evenly into the length of the window. Pad
        # with NaNs as needed to ensure it won't be skipped by the rolling range calculation below.
        # Take modulo to determine how much need to fill
        mod_incr = len(tr.data) % n_incr
        if mod_incr == 0:
            tr.data = np.append(tr.data, np.repeat(np.nan, n_samps - n_incr))
        else:
            tr.data = np.append(tr.data, np.repeat(np.nan, n_samps - mod_incr))

        # Calculate the rolling range
        Rrange = roll_range(tr.data, int(n_samps), int(n_incr), align="left")
        # Calculate max range
        maxRange = np.nanmax(Rrange)

        # Create dictionary output object
        data = {
            "snclq": snclq,
            "start_time": starttime,
            "end_time": endtime,
            "max_range": maxRange,
            "metric_name": "maxRangeMetric",
        }
        # Append trace output to list, will create a list of dictionaries
        d.append(data)

    # If database defined, insert metric information
    if database_config is not None:
        database = Database(**database_config)
        database.insert_metric(d)

    return d
