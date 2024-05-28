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

__all__ = ["DCOffSetTimesMetric"]

import numpy as np
from pycheron.util.logger import Logger
from pycheron.util.masks import samples2time
from pycheron.db.sqllite_db import Database


def DCOffSetTimesMetric(
    st,
    windowSecs=1800,
    incrementSecs=None,
    threshold=0.9,
    generateMasks=False,
    masksByTime=True,
    logger=None,
    database_config=None,
):
    """
    Metric to determine DC offset times (ie., where a shift in the signal mean is detected)

    :param st: ObsPy Stream object
    :type st: obspy.core.stream.Stream
    :param windowSecs: Chunk size (secs) used in DCOffset calculations (DEFAULT = 1800s)
    :type windowSecs: int
    :param incrementSecs: Increment (secs) for starttime of sequential chunks (DEFAULT = windowSecs/2)
    :type incrementSecs: int
    :param threshold: Threshold used for outlier detection
    :type threshold: float
    :param generateMasks: If True, a boolean mask will be created for dc offsets found.
    :type generateMasks: bool
    :param masksByTime: Boolean to determine whether masks are generated by time. If True, masks will be generated with
                        a start/end time, if false, they will be generated as boolean array.
    :type masksByTime: bool
    :param logger: (logger object) - If using a logger, (you must create one using the util.logger class)
    :type logger: pycheron.util.logger.Logger
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database

    :return: list of dictionaries (one for each trace in stream object) containing the following keys and types:

                * snclq (`str`): station, network, channel, location, quality information for each trace
                * start_time (`str`): start time of trace
                * end_time (`str`): end time of trace
                * dc_offset_times (`list`): times were dc offsets were found in the data
                * masks (`list` or `numpy.ndarray`) list of dictionaries of start/end times where dc offsets were found
                                                    or boolean array of 1's and 0's, where 1's indicate a dc offset
                * metric_name (`str`): metric name, DCOffsetTimesMetric

    :rtype: dict

    * Code originally ported from IRISMustangMetrics R Cran Package
    (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, March 20). CRAN-Package IRISMustangMetrics.
    The Comprehensive R Archive Network. Retrieved from
    https://cran.r-project.org/web/packages/IRISMustangMetrics/index.html) and augmented and adapted for use within
    Pycheron

    If difference in means between sequential chunks of signal is greater than the typical
    std deviation of a window then a DC offset has occurred
    (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)

    **Algorithm steps:**
    Steps taken directly from IRISMustangMetric documentation:
    (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)
    #. Break up signal into windows (dictated by windowSecs) spaced incrementSecs apart
    #. For each window calculate the signal mean and signal std deviation
    #. The resulting mean and std dev arrays are of length 47 for 24 hours of signal
    #. Metric = abs(lagged difference of window means)/mean(window std devs)
    #. DC offset = times when metric > threshold

    **Examples**

    .. code-block:: python

        # Read in data
        data = 'test/test_data/7A_CABN_ALL.988887.tar.mseed'
        st = obspy.read(data)
        # Calculate DCOffSetTimesmetric
        d = DCOffSetTimesMetric(st,1800,None,0.9, True)
        # Convert BHE stream to trace object
        tr = st[1]
        print 'dcOffsetTimes for BHE channel:', d[1]
        >>> 'dcOffsetTimes for BHE channel:'
            {'metric_name': 'DCOffsetTimesMetric',
            'dc_offset_times': [
                UTCDateTime(2013, 11, 1, 0, 45, 0, 25000),
                UTCDateTime(2013, 11, 1, 1, 15, 0, 25000),
                UTCDateTime(2013, 11, 1, 9, 30, 0, 25000),
                UTCDateTime(2013, 11, 1, 9, 45, 0, 25000),
                UTCDateTime(2013, 11, 1, 10, 0, 0, 25000),
                UTCDateTime(2013, 11, 1, 10, 15, 0, 25000),
                UTCDateTime(2013, 11, 1, 10, 30, 0, 25000),
                UTCDateTime(2013, 11, 1, 10, 45, 0, 25000),
                UTCDateTime(2013, 11, 1, 15, 45, 0, 25000),
                UTCDateTime(2013, 11, 1, 16, 15, 0, 25000)],
                'snclq': u'7A.CABN..BHE.M', 'start_time': UTCDateTime(2013, 11, 1, 0, 0),
                'masks': array([0, 0, 1, ..., 0, 0, 0]), 'end_time': UTCDateTime(2013, 11, 1, 23, 59, 59, 975000)}

    .. code-block:: python

        # Plotting
        de = []
        for i in d[1]['dc_offset_times']:
            de.append(timedelta.total_seconds(i.datetime-tr.stats.starttime.datetime))
        plt.plot(tr.times(),tr.data)
        for i in de:
            plt.axvline(x=i,color='black')
        plt.xlabel('Seconds from %s' %str(tr.stats.starttime))
        plt.xlim([min(tr.times()),max(tr.times())])
        plt.ylabel('Amplitude')
        plt.title('%s WindowSecs 1800, incrementSecs 900, threshold 0.9'%str(tr.id))
        plt.show()

    .. image:: _static/DCOffSetTimesMetric.png

    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Initialize incrementSecs variable to default value if None
    if incrementSecs is None:
        incrementSecs = windowSecs / 2

    # Initialize d, a list of dictionaries to fill with metric values
    d = []

    # Merge stream objects with same IDs, if gaps present, fill with '999999', then remove '999999' values from trace
    # Make copy of stream object to retain original stream, otherwise will be permanently changed
    stC = st.copy()
    stM = stC.merge(fill_value=999999)
    # May want to rethink this step of removing fill values
    for trace in stM:
        tr = trace
        tr.data = tr.data[tr.data != 999999]

    # Loop through traces and get trace information, window information, then calculate mean and std dev for each window
    for trace in stM:
        # Get trace information: sample rate, snclq, starttime, endtime for each trace
        tr = trace
        fs = tr.stats.sampling_rate

        # Get sncl
        snclq = tr.get_id()

        # Get starttime, endtime
        starttime = tr.stats.starttime
        endtime = tr.stats.endtime

        # Calculate window information
        windowSamples = windowSecs * tr.stats.sampling_rate
        incrementSamples = incrementSecs * tr.stats.sampling_rate
        outLength = int(len(tr) / incrementSamples)

        # Initialize mean, standard deviation, and indices lists
        means = []
        sds = []
        indices = []

        # Loop through all of the window chunks and calculate the mean and std dev
        for i in range(outLength):
            lo = int(i * incrementSamples + 1)
            hi = int(lo + windowSamples)
            means.append(np.mean(tr.data[lo:hi]))
            sds.append(np.std(tr.data[lo:hi]))
            indices.append(int(lo))

        # Calculate the difference between sequential means
        # Determine indices where the value is greater than the provided threshold
        # Use those indices to determine dcOffset times and append to dcOffsetTimes list
        metric = abs(np.diff(means, n=1, axis=0)) / np.mean(sds)
        jumps = [index for index, value in enumerate(metric) if value > threshold]
        # Initialize dcOffsetTimes list to append list of dcOffsettimes to
        dcOffsetTimes = []
        for i in jumps:
            times = starttime + indices[i + 1] / tr.stats.sampling_rate
            dcOffsetTimes.append(times.isoformat())

        # Create masks is generateMasks = True, first initialize masks to None
        masks = None
        if generateMasks and jumps:
            # If masksByTime True, then get masks by time
            if masksByTime:
                masks = samples2time(np.asarray(jumps), fs, starttime)
            # Otherwise do boolean
            else:
                masks = np.zeros(len(tr.data), dtype=int)
                masks[jumps] = 1

        # Create metrics list to append to from each stream object
        metrics = {
            "snclq": snclq,
            "start_time": starttime.isoformat(),
            "end_time": endtime.isoformat(),
            "dc_offset_times": dcOffsetTimes,
            "masks": masks,
            "metric_name": "DCOffsetTimesMetric",
        }

        d.append(metrics)

    # If database defined, insert metric information
    if database_config is not None:
        database = Database(**database_config)
        database.insert_metric(d)

    return d
