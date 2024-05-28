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
from pycheron.util.masks import consecutive
from pycheron.db.sqllite_db import Database

__all__ = ["deadChanMean"]


def deadChanMean(
    st,
    win_size=None,
    threshold=0.05,
    generateMasks=False,
    masksByTime=True,
    database_config=None,
):
    """
    Function uses a simplistic rolling mean, using win_size and increments over the trace object by that win_size.
    Unlike a typical rolling mean, no samples are shared between windows, e.g.

    `[-----------][-----------][-----------][-----------][-----------][-----------][-----------][-----------]`

    :param st: - ObsPy stream object
    :type st: obspy.core.stream.Stream
    :param win_size: number of samples to use for windowing. When None, determines the number of samples to use based on
                     the sample rate
    :type win_size: int
    :param threshold: threshold to trigger deadChan for within rolling window. Recommended to keep this default value
                      for triggering a dead channel
    :type threshold: int
    :param generateMasks: Generate boolean qc mask array.
    :type generateMasks: bool
    :param masksByTime: Generate a time-based qc mask list.
    :type masksByTime: bool
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database

    :return: Returns list of dictionaries that includes the following keys and type for each trace within the stream:

            * dead_channel_times (`list`) - Lists of dicts with dead channel start/endtime
            * end_time (`str`) - end time of trace
            * start_time (`str`) - start time of trace
            * masks (`list` or `numpy.ndarray`) - masks for that trace, which is either a boolean array or a list of
                                                start/stop times in UTC
            * metric_name (`str`) - deadChanADFMetric string
            * snlcq (`str`) - station network location channel information for the specific trace

    :rtype: list



    """

    # Loop through stream object and calculate a moving average for each trace based on user defined window size and
    # period -- Using pandas built in moving average because very fast

    out = []
    for s in range(len(st)):
        # Pull out trace one by one and loop through
        tr = st[s]

        # Initialize list to hold dad channel times
        dead_times = []
        # Obtain appropriate window size based on the sample rate
        if win_size is None:
            if int(tr.stats.sampling_rate) == 100:
                win_size = 28000
            elif int(tr.stats.sampling_rate) == 40:
                win_size = 11200
            elif int(tr.stats.sampling_rate) == 20:
                win_size = 5600
            # TODO: May want to think about changing this size for other sample rates -- maybe split into higher than
            # TODO: 100sps and lower than 20sps?
            else:
                win_size = 30000

        # Initialize empty list to append means to
        rolling_mean = []

        # Loop through data and calculate the mean of each window
        for n in np.arange(0, len(tr.data) - win_size, win_size):
            # Define n_start, n_end
            n_start = n
            n_end = n + win_size

            # Calculate rolling_mean
            rolling_mean.append(np.mean(tr.data[n_start:n_end]))

        # Convert list to array
        rolling_mean = np.asarray(rolling_mean)

        # Get differences in consecutive array elements
        diffs = np.ediff1d(rolling_mean)

        # Find indices where differences are greater than threshold
        inds = np.argwhere(abs(diffs) <= threshold).ravel()

        # if inds returns an empty list, skip everything
        if inds.size != 0:
            # Convert indices to time for masks.
            # If consecutive indices, use start/end time. If only one use ind as start
            # and end time. E.g., for long periods of time use first index and last index until a break in being below
            # the threshold

            # Group array into consecutive element arrays to identify stop/start times easier
            cons = consecutive(inds)

            # Loop through range of consecutive elements to get the start/end time of dead Channel issues
            for c in range(len(cons)):
                # Skip single point elements as these are likely flukes and only process consecutive elements
                if len(cons[c]) != 1:
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
                    # Loop through consecutive element blocks, skip single element issues. Append Starttime and Endtime
                    # to list, then convert to array
                    for c in range(len(cons)):
                        if len(cons[c]) != 1:
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
                        "metric_name": "deadChanMeanMetric",
                    }
                    out.append(dt)
                else:
                    # if maskByTime is false create binary array
                    deadMask = np.zeros(len(tr.data), dtype=int)
                    for c in range(len(cons)):
                        if len(cons[c]) != 1:
                            deadMask[cons[c] * win_size] = 1
                    dt = {
                        "snclq": tr.get_id(),
                        "start_time": tr.stats.starttime.isoformat(),
                        "end_time": tr.stats.endtime.isoformat(),
                        "dead_channel_times": dead_times if dead_times else None,
                        "masks": deadMask,
                        "metric_name": "deadChanMeanMetric",
                    }
                    out.append(dt)

            else:

                dt = {
                    "snclq": tr.get_id(),
                    "start_time": tr.stats.starttime.isoformat(),
                    "end_time": tr.stats.endtime.isoformat(),
                    "dead_channel_times": dead_times if dead_times else None,
                    "masks": "No Masks Created",
                    "metric_name": "deadChanMeanMetric",
                }
                out.append(dt)
        else:
            dt = {
                "snclq": tr.get_id(),
                "start_time": tr.stats.starttime.isoformat(),
                "end_time": tr.stats.endtime.isoformat(),
                "dead_channel_times": None,
                "masks": None,
                "metric_name": "deadChanMeanMetric",
            }
            out.append(dt)

    # If database defined, insert metric information
    if database_config is not None:
        database = Database(**database_config)
        database.insert_metric(out)

    return out
