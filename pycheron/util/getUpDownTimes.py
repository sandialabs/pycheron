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

__all__ = ["getUpDownTimes"]

from pycheron.util.getUniqueIds import getUniqueIds
from pycheron.util.logger import Logger


def getUpDownTimes(st, min_signal=180, min_gap=300, logger=None):
    """

    Obtain on/off times of a ObsPy Stream object. Traces with a duration less than min_signal are ignored. Gaps less
    than min_gap are also ignored. On times will be continuous segments, while off-times will represent gaps in the
    data.

    :param st: ObsPy Stream object
    :type st: obspy.core.stream.Stream
    :param min_signal: Minimum trace duration in seconds.
    :type min_signal: int
    :param min_gap: Minimum gap in seconds.
    :type min_gap: int

    :return: Returns a list of dictionaries with the following keys and types:

            * id (`str`)
            * start_time (`list`)
            * end_time (`list`)
            * channel_up_time (`list`)

    :rtype: list

    * Code originally ported from IRISSeismic R Cran Package
      (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
      The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index) and
      augmented and adapted for use within Pycheron.

    **Example**

    .. code-block:: python

        import obspy
        from pycheron.util.getUpDownTimes import getUpDownTimes

        # test data
        data = 'test/test_data/6e_sp06_ehe.407438.tar.mseed'

        # reading in stream
        st = obspy.read(data)


        # Defining parameters-- these are optional, defaults are: min_sig=30; min_gaps=60
        min_sig = 30
        min_gaps = 60

        # getting up down times, returns a list with dictionaries in it. If there is only 1 trace it seems redundant,
        # but the function is made to handle streams with multiple traces.
        upDown_one = getUpDownTimes(st_one,min_sig,min_gaps)
        >>> [{'Channel Up-time': 0, 'Endtime': UTCDateTime(2014, 1, 18, 23, 59, 59, 995000), 'id': u'6E.SP06..EHE',
              'Starttime': UTCDateTime(2014, 1, 1, 0, 0, 0, 5000)}]


    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # getting list of unique identifiers
    u_ids = getUniqueIds(st)

    # Initialize lists
    up_down_times_list = []
    starttimes = []
    endtimes = []

    # If there are many unique identifiers
    if len(u_ids) > 1:

        # looping through trace in stream
        for trace in st:

            # Get start/end times
            starttime = trace.stats.starttime
            endtime = trace.stats.endtime

            # determine signal duration for each trace
            signal_duration = endtime - starttime
            if signal_duration >= min_signal:
                starttimes.append(starttime)
                endtimes.append(endtime)
                duration = signal_duration
            else:
                logger.warn("getUpDownTimes warning: Signal Duration < Minimum Signal")
                continue

            # Fill out dictionary with id, start time, end time, and channel up time
            up_down_times = {
                "id": trace.get_id(),
                "start_time": starttimes,
                "end_time": endtimes,
                "channel_up_time": duration,
            }
            # get the difference in seconds between the current end time
            # and the next start time.

            # If there are only 2 traces in stream
            if len(starttimes) == 2:
                delta = endtimes[0] - starttimes[1]

                # if data is below threshold, continuous
                if delta < 0 or delta < min_gap:
                    up_down_times["start_time"] = "NaN"
                    up_down_times["end_time"] = "NaN"
                    up_down_times["channel_up_time"] = 100
                    up_down_times_list.append(up_down_times)

                # If data is above threshold, discontinuous
                else:
                    up_down_times["start_time"] = endtimes[0]
                    up_down_times["end_time"] = starttimes[1]
                    up_down_times["channel_up_time"] = delta
                    up_down_times_list.append(up_down_times)

            # If there are more than two traces
            else:
                # looping through traces
                for i in range(len(starttimes) - 1):
                    delta = endtimes[i] - starttimes[i + 1]
                    # if data is below threshold, continuous
                    if delta < 0 or delta < min_gap:
                        up_down_times["start_time"] = "NaN"
                        up_down_times["end_time"] = "NaN"
                        up_down_times["channel_up_time"] = 100
                        up_down_times_list.append(up_down_times)

                    # If data is above threshold, discontinuous
                    else:
                        up_down_times["start_time"] = endtimes[i]
                        up_down_times["end_time"] = starttimes[i + 1]
                        up_down_times["channel_up_time"] = delta
                        up_down_times_list.append(up_down_times)

        return up_down_times_list

    # If there is only one unique identifier
    elif len(u_ids) == 1:

        # Create empty list for start/end times
        starttimes = []
        endtimes = []
        # looping through unique stream
        for trace in st:
            # Get start/end times for each trace
            starttime = trace.stats.starttime
            endtime = trace.stats.endtime

            # determine signal duration for each trace
            signal_duration = endtime - starttime
            if signal_duration >= min_signal:
                starttimes.append(starttime)
                endtimes.append(endtime)
                duration = signal_duration
            else:
                print("getUpDownTimes Error: Signal Duration < Minimum Signal")
                logger.log("getUpDownTimes Error: Signal Duration < Minimum Signal")
                continue

            # Fill out dictionary with id, start time, end time, and channel up time
            up_down_times = {
                "id": trace.get_id(),
                "start_time": starttime,
                "end_time": endtime,
                "channel_up_time": duration,
            }
            # get the difference in seconds between the current end time
            # and the next start time.
            if len(starttimes) == 1:
                up_down_times_list.append(up_down_times)
                up_down_times["channel_up_time"] = 100  # was 0

            # If there is only 2 traces in stream
            elif len(starttimes) == 2:
                delta = endtimes[0] - starttimes[1]

                # if data is below threshold, continuous
                if delta < 0 or delta < min_gap:
                    up_down_times["start_time"] = "NaN"
                    up_down_times["end_time"] = "NaN"
                    up_down_times["channel_up_time"] = 100
                    up_down_times_list.append(up_down_times)

                # If data is above threshold, discontinuous
                else:
                    up_down_times["start_time"] = endtimes[0]
                    up_down_times["end_time"] = starttimes[1]
                    up_down_times["channel_up_time"] = delta
                    up_down_times_list.append(up_down_times)

            # If there are more than two traces
            else:
                # looping through traces
                for i in range(len(starttimes) - 1):
                    delta = endtimes[i] - starttimes[i + 1]
                    # if data is below threshold, continuous
                    if delta < 0 or delta < min_gap:
                        up_down_times["start_time"] = "NaN"
                        up_down_times["end_time"] = "NaN"
                        up_down_times["channel_up_time"] = 100
                        up_down_times_list.append(up_down_times)

                        # If data is above threshold, discontinuous
                    else:
                        up_down_times["start_time"] = endtimes[i]
                        up_down_times["end_time"] = starttimes[i + 1]
                        up_down_times["channel_up_time"] = delta
                        up_down_times_list.append(up_down_times)

        return up_down_times_list
