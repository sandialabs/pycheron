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

__all__ = ["gapMetric"]

import numpy as np
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.trace import Trace
from pycheron.util.masks import samples2time
from pycheron.util.logger import Logger


def gapMetric(
    st,
    generateMasks=False,
    separateMasks=True,
    completeDay=True,
    masksByTime=True,
    logger=None,
    database=None,
):
    """
    Function calculates metrics associated with gaps and overlaps in a seismic signal.

    :param st: obspy stream object
    :type st: obspy.core.stream.Stream
    :param generateMasks: (boolean) If True, QC masks will be created.
    :type generateMasks: bool
    :param separateMasks: (boolean) If True, 2 masks will be create. 1 for gaps, and 1 for overlaps, if False these masks will be merged and 3 masks will be returned. Gaps, Overlaps, Combined.
    :type separateMasks: bool
    :param completeDay: (boolean) If True, if an individual day is does not start at 00:00:00 or end at 23:59:59, a small trace will be made to complete the day. This should be used if, for example, have only 6 hours of data, but would like to consider that 6 hours out of a whole day, thus gaps
                               would be found and masks would be created around that segment of data.
    :type completeDay: bool
    :param masksByTime: (bool) - Boolean to determine whether masks are generated by time. If True, masks will be generated with a start/end time, if false, they will be generated as boolean array.
    :type masksByTime: bool
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database

    :return: Tuple of station gap summary information and channel gap information. See below.
    :rtype: tuple

    .. code-block:: python

        # Example, returns 2 objects. Each object is a list of dictionaries
        summary, data = gapMetric()

    **Metrics include:**

    * Station Summary data:
        * `station_completeness` - total percentage of total requested time for which signal is available for a station (average of channel percentages)
        * `network` - network name
        * `station` - station name
        * `location` - location
        * `total_gaps` - total number of gaps in stream
        * `total_overlaps` - total number of overlaps in stream
        * `metric_name` - gapMetricStation
        * `starttime` - end time of trace
        * `endtime` - start time of trace
        * `gap_masks` - Mask for gapped data. Boolean array
        * `overlap_masks` - Mask for overlapped data. Boolean array
        * `combined_masks` - Mask for gapped and overlaped data. Boolean array
        * `maximum_gap` - length of maximum gap found in channel within stream
        * `maximum_overlap` - length of maximum overlap(s) found in channel within stream
        * `channel_percent_available` - percentage  of total requested time for which a signal is available for channel

    * Channel data
        * `channel` - channel that the above metrics apply to
        * `duration` - duration of gaps/overlaps within channel
        * `type` - overlap or gap
        * `gap_overlap_start_time` - gap/overlap start time
        * `gap_overlap_end_time` - gap/overlap end time
        * `metric_name` - gapChannelMetric

    **Examples**

    1. Gaps with multiple channels, same station

    .. code-block:: python

        #test data with gaps with multiple channels, same station
        data = 'test/test_data/6e_sp06_Eall.397679.tar.mseed'
        #reading in stream
        st = obspy.read(data)
        gaps, sum = gapMetric(st)
        print 'Gap Details:'
        for i in gaps:
            print i
        >>> Gap Details:
        >>> [{'duration': 408.0050001144409, 'gap/overlap_start_time': '2014-01-19T21:01:07.985000', 'gap/overlap_end_time': '2014-01-19T21:07:56', 'type': 'Gap', 'channel': u'EHE'}, {'duration': 86397.97009897232, 'gap/overlap_start_time': '2014-01-20T00:00:00.990000', 'gap/overlap_end_time': '2014-01-20T23:59:58.970099', 'type': 'Gap', 'channel': u'EHE'}]
        >>> [{'duration': 408.0050001144409, 'gap/overlap_start_time': '2014-01-19T21:01:07.985000', 'gap/overlap_end_time': '2014-01-19T21:07:56', 'type': 'Gap', 'channel': u'EHN'}, {'duration': 86397.97009897232, 'gap/overlap_start_time': '2014-01-20T00:00:00.990000', 'gap/overlap_end_time': '2014-01-20T23:59:58.970099', 'type': 'Gap', 'channel': u'EHN'}]
        >>> [{'duration': 408.0050001144409, 'gap/overlap_start_time': '2014-01-19T21:01:07.985000', 'gap/overlap_end_time': '2014-01-19T21:07:56', 'type': 'Gap', 'channel': u'EHZ'}, {'duration': 86397.97009897232, 'gap/overlap_start_time': '2014-01-20T00:00:00.990000', 'gap/overlap_end_time': '2014-01-20T23:59:58.970099', 'type': 'Gap', 'channel': u'EHZ'}]

        print "Gap Summary:"
        for i in sum:
            print i
        >>> Gap Summary:
        >>> {'overlap_masks': None, 'station_completeness': 49.76476708626203, 'total_gaps': 6, 'start_time': UTCDateTime(2014, 1, 19, 0, 0, 0, 5000), 'total_overlaps': 0, 'combined_masks': None, 'maximum_overlap': nan, 'channel_percent_available': 49.76476708626204, 'end_time': UTCDateTime(2014, 1, 20, 23, 59, 58, 995000), 'gap_masks': None, 'maximum_gap': 86397.97009897232, 'channel': u'EHE', 'metirc_name': 'gapMetric'}
        >>> {'overlap_masks': None, 'station_completeness': 49.76476708626203, 'total_gaps': 6, 'start_time': UTCDateTime(2014, 1, 19, 0, 0, 0, 5000), 'total_overlaps': 0, 'combined_masks': None, 'maximum_overlap': nan, 'channel_percent_available': 49.76476708626204, 'end_time': UTCDateTime(2014, 1, 20, 23, 59, 58, 995000), 'gap_masks': None, 'maximum_gap': 86397.97009897232, 'channel': u'EHN', 'metirc_name': 'gapMetric'}
        >>> {'overlap_masks': None, 'station_completeness': 49.76476708626203, 'total_gaps': 6, 'start_time': UTCDateTime(2014, 1, 19, 0, 0, 0, 5000), 'total_overlaps': 0, 'combined_masks': None, 'maximum_overlap': nan, 'channel_percent_available': 49.76476708626204, 'end_time': UTCDateTime(2014, 1, 20, 23, 59, 58, 995000), 'gap_masks': None, 'maximum_gap': 86397.97009897232, 'channel': u'EHZ', 'metirc_name': 'gapMetric'}

    2. Data with overlaps single channel

    .. code-block:: python

        #test data with overlaps single channel
        data = 'test/test_data/qualityflags.mseed'
        st = obspy.read(data)
        gaps, sum= gapMetric(st)
        print 'Gap Details:'
        >>> Gap Details:
        >>> [{'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': -2.054999828338623, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2007-12-31T23:59:59.915000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': 86397.01009917259, 'gap/overlap_start_time': '2008-01-01T00:00:01.970000', 'gap/overlap_end_time': '2008-01-01T23:59:58.985099', 'type': 'Gap', 'channel': u'EHE'}]

        for i in gaps:
            print i
        print "Gap Summary:"
        for i in sum:
            print i
        >>> Gap Summary:
        >>> {'overlap_masks': None, 'station_completeness': 0.0023957440604789326, 'total_gaps': 1, 'start_time': UTCDateTime(2007, 12, 31, 23, 59, 59, 915000), 'total_overlaps': 17, 'combined_masks': None, 'maximum_overlap': 2.054999828338623, 'channel_percent_available': 0.0023957440604789326, 'end_time': UTCDateTime(2008, 1, 1, 23, 59, 58, 995000), 'gap_masks': None, 'maximum_gap': 86397.01009917259, 'channel': u'EHE', 'metirc_name': 'gapMetric'}

    3. Data with gaps and overlaps single channel

    .. code-block:: python

        #test data with gaps and overlaps single channel
        data = 'test/test_data/sglchan_onegap_oneoverlap.mseed'
        st = obspy.read(data)
        gaps, sum= gapMetric(st)
        print 'Gap Details:'
        for i in gaps:
            print i
        >>> Gap Details:
        >>> [{'duration': 16.479999780654907, 'gap/overlap_start_time': '2008-01-01T00:01:34.480000', 'gap/overlap_end_time': '2008-01-01T00:01:50.965000', 'type': 'Gap', 'channel': u'EHE'}, {'duration': -2.06000018119812, 'gap/overlap_start_time': '2008-01-01T00:01:53.020000', 'gap/overlap_end_time': '2008-01-01T00:01:50.965000', 'type': 'Overlap', 'channel': u'EHE'}, {'duration': 86191.2000989914, 'gap/overlap_start_time': '2008-01-01T00:03:27.780000', 'gap/overlap_end_time': '2008-01-01T23:59:58.985099', 'type': 'Gap', 'channel': u'EHE'}]

        print "Gap Summary:"
        for i in sum:
            print i
        >>> Gap Summary:
        >>> {'overlap_masks': None, 'station_completeness': 0.22170325039695626, 'total_gaps': 2, 'start_time': UTCDateTime(2007, 12, 31, 23, 59, 59, 765000), 'total_overlaps': 1, 'combined_masks': None, 'maximum_overlap': 2.06000018119812, 'channel_percent_available': 0.22170325039695626, 'end_time': UTCDateTime(2008, 1, 1, 23, 59, 58, 995000), 'gap_masks': None, 'maximum_gap': 86191.2000989914, 'channel': u'EHE', 'metirc_name': 'gapMetric'}

    """
    if logger == None:
        logger = Logger(None)

    # ------------------Initialize---------------------------------------------
    # Initalize list of counts of various metrics on gaps/overlaps
    max_gap_listE = []
    max_overlap_listE = []
    max_gap_listN = []
    max_overlap_listN = []
    max_gap_listZ = []
    max_overlap_listZ = []

    # Initalize channel lists, totalSec, duration, and cont_chunk
    chans_list = []
    uniq_chans = []
    totalSec = []
    durSeg = []

    # Initialize output variable data
    data = []

    # Initialize percent available list
    percent_avail_list = []
    sta_percent = []
    percent_av = None

    # Initialize values
    num_gaps = 0
    max_gap = 0
    num_overlaps = 0
    max_overlap = 0
    num_gaps_per_sta = 0
    num_overlaps_per_sta = 0
    total_gaps = 0
    total_overlaps = 0
    sampRate = 0
    gapMask = None
    overlapMask = None
    allMasks = None
    gStarttime = []
    gEndtime = []
    num_segs = []

    # Loop through traces get channel list (will be repeat) and number of segments

    trM = st.copy()

    if completeDay == True:
        trM = trM.merge()

        for i in range(len(trM)):
            tr = trM[i]
            bod = UTCDateTime(
                tr.stats.endtime.year,
                tr.stats.endtime.month,
                tr.stats.endtime.day,
                00,
                00,
                00,
                00,
            )
            eod = UTCDateTime(
                tr.stats.endtime.year,
                tr.stats.endtime.month,
                tr.stats.endtime.day,
                23,
                59,
                59,
                99,
            )
            end = tr.stats.endtime
            start = tr.stats.starttime
            if end <= eod:
                tr_data = np.repeat(np.array([99999], dtype=np.int32), 3)
                header = {
                    "sampling_rate": tr.stats.sampling_rate,
                    "calib": tr.stats.calib,
                    "npts": 3,
                    "network": tr.stats.network,
                    "station": tr.stats.station,
                    "channel": tr.stats.channel,
                    "location": tr.stats.location,
                    "starttime": eod - (3.0 / tr.stats.sampling_rate),
                    "endtime": eod,
                }
                trNew = Trace(tr_data, header)
                trM = trM.append(trNew)
                logger.log("gapMetric(): New Trace Added")

            if start > bod:
                tr_data = np.repeat(np.array([99999], dtype=np.int32), 3)
                header = {
                    "sampling_rate": tr.stats.sampling_rate,
                    "calib": tr.stats.calib,
                    "npts": 3,
                    "network": tr.stats.network,
                    "station": tr.stats.station,
                    "channel": tr.stats.channel,
                    "location": tr.stats.location,
                    "starttime": bod,
                    "endtime": bod + (3.0 / tr.stats.sampling_rate),
                }
                trNew = Trace(tr_data, header)
                trM = trM.append(trNew)
                logger.log("gapMetric(): New Trace Added")

    else:
        # Get total time duration for each stream object by merging all gaps together
        trM = trM.merge()

    for i in range(len(st)):
        tr = st[i]
        network = tr.stats.network
        station = tr.stats.station
        location = tr.stats.location
        chans_list.append(tr.stats.channel)
        # Get number of segments per channel
    num_segs = [chans_list.count(chan) for chan in np.unique(chans_list)]

    fs = st[0].stats.sampling_rate

    for i in range(len(trM)):
        tr = trM[i]
        totalSec.append(tr.stats.endtime - tr.stats.starttime)
        uniq_chans.append(tr.stats.channel)

    # For multiple gaps per station create list of totSec to iterate through that are the same number of iterations as
    # gaps, also make one for channel list
    totSec = []
    chans = []
    for i in range(len(num_segs)):
        totSec.append(np.repeat(totalSec[i], num_segs[i] - 1))
        chans.append(np.repeat(uniq_chans[i], num_segs[i] - 1))
    totSec = np.concatenate(totSec)
    chans = np.concatenate(chans)
    chans.sort()

    # ------------------Gap Function-------------------------------------------
    gaps = st.get_gaps()

    index_gaps = np.empty(0, dtype=int)
    index_over = np.empty(0, dtype=int)

    for i in range(len(gaps)):
        # ------------------------QC MASKS----------------------
        # create gap masks
        if gaps[i][7] > 0:
            start_gap = gaps[i][4]
            secFromStart = start_gap - trM[0].stats.starttime
            maskStartIndex = np.floor(secFromStart * sampRate)
            maskEndIndex = maskStartIndex + gaps[i][7]
            index = np.arange(maskStartIndex, maskEndIndex, 1, dtype=int)
            index_gaps = np.append(index_gaps, index)
        # create overlap masks
        if gaps[i][7] < 0:
            start_over = gaps[i][5]
            secFromStart = start_over - trM[0].stats.starttime
            maskStartIndex = np.abs(np.floor(secFromStart * sampRate))
            maskEndIndex = maskStartIndex + np.abs(gaps[i][7])
            index = np.arange(maskStartIndex, maskEndIndex, 1, dtype=int)
            index_over = np.append(index_over, index)
        # -------------------------------------------------------------

        if gaps[i][3].endswith("E") or gaps[i][3].endswith("2"):
            if gaps[i][6] > 0:
                max_gap_listE.append(gaps[i][6])
                max_gap = max(max_gap_listE)
                num_gaps_per_sta = len(max_gap_listE)
            elif gaps[i][6] < 0:
                max_overlap_listE.append(abs(gaps[i][6]))
                max_overlap = max(max_overlap_listE)
                num_overlaps_per_sta = len(max_overlap_listE)
                percent_avail_list.append(100)
            # getting gap/overlap start and endtime
            if num_gaps_per_sta != 0 or num_overlaps_per_sta != 0:
                gStarttime.append(gaps[i][4])
                gEndtime.append(gaps[i][5])
            else:
                gStarttime = None
                gEndtime = None
            durSeg.append(gaps[i][6])

        elif gaps[i][3].endswith("N") or gaps[i][3].endswith("1"):
            if gaps[i][6] > 0:
                max_gap_listN.append(gaps[i][6])
                max_gap = max(max_gap_listN)
                num_gaps_per_sta = len(max_gap_listN)
            elif gaps[i][6] < 0:
                max_overlap_listN.append(abs(gaps[i][6]))
                max_overlap = max(max_overlap_listN)
                num_overlaps_per_sta = len(max_overlap_listN)
                percent_avail_list.append(100)
            # getting gap/overlap start and endtime
            if num_gaps_per_sta != 0 or num_overlaps_per_sta != 0:
                gStarttime.append(gaps[i][4])
                gEndtime.append(gaps[i][5])
            else:
                gStarttime = None
                gEndtime = None
            durSeg.append(gaps[i][6])

        elif gaps[i][3].endswith("Z"):
            if gaps[i][6] > 0:
                max_gap_listZ.append(gaps[i][6])
                max_gap = max(max_gap_listZ)
                num_gaps_per_sta = len(max_gap_listZ)
            elif gaps[i][6] < 0:
                max_overlap_listZ.append(abs(gaps[i][6]))
                max_overlap = max(max_overlap_listZ)
                num_overlaps_per_sta = len(max_overlap_listZ)
                percent_avail_list.append(100)
            # getting gap/overlap start and endtime
            if num_gaps_per_sta != 0 or num_overlaps_per_sta != 0:
                gStarttime.append(gaps[i][4])
                gEndtime.append(gaps[i][5])
            else:
                gStarttime = None
                gEndtime = None
            durSeg.append(gaps[i][6])

        if i == len(gaps) - 1 or chans[i] != chans[i + 1]:
            ind_gaps = []
            if gaps[i][3].endswith("E") or gaps[i][3].endswith("2"):
                if len(gaps) == 1:
                    if gaps[i][6] > 0:
                        percent_av = 100 - (100 * max_gap_listE[0] / totalSec[0])
                    else:
                        percent_av = 100
                else:
                    gap_secsE = sum(max_gap_listE)
                    if len(uniq_chans) == 1:
                        percent_avail_list.append(100 - (100 * gap_secsE / totSec[0]))
                    else:
                        percent_avail_list.append(100 - (100 * gap_secsE / totSec[i]))
                    percent_av = [x for x in percent_avail_list if x != 100]
                    if not percent_av:
                        percent_av = 100

            elif gaps[i][3].endswith("N") or gaps[i][3].endswith("1"):
                if len(gaps) == 1:
                    if gaps[i][6] > 0:
                        percent_av = 100 - (100 * max_gap_listN[0] / totalSec[0])
                    else:
                        percent_av = 100
                else:
                    gap_secsN = sum(max_gap_listN)
                    if len(uniq_chans) == 1:
                        percent_avail_list.append(100 - (100 * gap_secsN / totSec[0]))
                    else:
                        percent_avail_list.append(100 - (100 * gap_secsN / totSec[i]))
                    percent_av = [x for x in percent_avail_list if x != 100]
                    if not percent_av:
                        percent_av = 100

            elif gaps[i][3].endswith("Z"):
                if len(gaps) == 1:
                    durSeg = gaps[0][6]
                    if gaps[i][6] > 0:
                        percent_av = 100 - (100 * max_gap_listZ[0] / totalSec[0])
                    else:
                        percent_av = 100
                else:
                    gap_secsZ = sum(max_gap_listZ)
                    if len(uniq_chans) == 1:
                        percent_avail_list.append(100 - (100 * gap_secsZ / totSec[0]))
                    else:
                        percent_avail_list.append(100 - (100 * gap_secsZ / totSec[i]))
                    percent_av = [x for x in percent_avail_list if x != 100]
                    if not percent_av:
                        percent_av = 100
            try:
                for j in range(len(durSeg)):
                    d = {
                        "channel": chans[i],
                        "duration": durSeg[j],
                        "gap_overlap_start_time": gStarttime[j].isoformat(),
                        "gap_overlap_end_time": gEndtime[j].isoformat(),
                        "type": "Overlap" if durSeg[j] < 0 else "Gap",
                        "metric_name": "gapMetricChannel",
                    }
                    ind_gaps.append(d)
            except TypeError:
                d = {
                    "channel": chans[i],
                    "duration": durSeg,
                    "gap_overlap_start_time": gStarttime[0].isoformat(),
                    "gap_overlap_end_time": gEndtime[0].isoformat(),
                    "type": "Overlap" if durSeg < 0 else "Gap",
                    "metric_name": "gapMetricChannel",
                }
                ind_gaps.append(d)
            if ind_gaps:
                data.append(ind_gaps)

            sta_percent.append(percent_av)
            total_gaps += num_gaps_per_sta
            total_overlaps += num_overlaps_per_sta
            durSeg = []

        # station summary
        if generateMasks == True:
            if not masksByTime:
                index_gaps = np.unique(index_gaps)
                index_over = np.unique(index_over)
                gapMask = np.zeros(len(trM[0].data), dtype=int)
                gapMask[index_gaps] = 1

                overlapMask = np.zeros(len(trM[0].data), dtype=int)
                overlapMask[index_over] = 1
                if separateMasks == False:
                    allMasks = np.zeros(len(trM[0].data), dtype=int)
                    allMasks[index_gaps] = 1
                    allMasks[index_over] = 1

            else:
                index_gaps = np.unique(index_gaps)
                index_over = np.unique(index_over)

                gapMask = samples2time(index_gaps, fs, trM[0].stats.starttime)
                overlapMask = samples2time(index_over, fs, trM[0].stats.starttime)

                if separateMasks == False:
                    allMasks = samples2time(index_gaps, fs, trM[0].stats.starttime)
                    allMasks = allMasks + samples2time(
                        index_over, fs, trM[0].stats.starttime
                    )

    sum_data = []
    if len(data) == 0:
        additions = {
            "station_completeness": 100,
            "network": network,
            "station": station,
            "location": location,
            "total_gaps": total_gaps,
            "total_overlaps": num_overlaps,
            "metric_name": "gapMetricStation",
            "start_time": trM[0].stats.starttime,
            "end_time": trM[0].stats.endtime,
            "gap_masks": gapMask,
            "overlap_masks": overlapMask,
            "combined_masks": allMasks,
            "channel_percent_available": 100,
        }
        sum_data.append(additions)
    else:
        for k in range(len(data)):
            if data[k][0]["channel"].endswith("E") or data[k][0]["channel"].endswith(
                "1"
            ):
                if max_overlap_listE:
                    max_overlap = np.max(max_overlap_listE)
                else:
                    max_overlap = np.nan
                if max_gap_listE:
                    max_gap = np.max(max_gap_listE)
                else:
                    max_gap = np.nan
                try:
                    channel_per = percent_av[k]
                except TypeError:
                    channel_per = percent_av

            elif data[k][0]["channel"].endswith("N") or data[k][0]["channel"].endswith(
                "2"
            ):
                if max_overlap_listE:
                    max_overlap = np.max(max_overlap_listN)
                else:
                    max_overlap = np.nan
                if max_gap_listE:
                    max_gap = np.max(max_gap_listN)
                else:
                    max_gap = np.nan
                try:
                    channel_per = percent_av[k]
                except TypeError:
                    channel_per = percent_av
            elif data[k][0]["channel"].endswith("Z"):
                if max_overlap_listE:
                    max_overlap = np.max(max_overlap_listZ)
                else:
                    max_overlap = np.nan
                if max_gap_listE:
                    max_gap = np.max(max_gap_listZ)
                else:
                    max_gap = np.nan
                try:
                    channel_per = percent_av[k]
                except TypeError:
                    channel_per = percent_av

            if len(data[k]) == 0:
                additions = {
                    "station_completeness": 100,
                    "network": network,
                    "station": station,
                    "location": location,
                    "total_gaps": total_gaps,
                    "total_overlaps": num_overlaps,
                    "metric_name": "gapMetricStation",
                    "start_time": trM[0].stats.starttime,
                    "end_time": trM[0].stats.endtime,
                    "gap_masks": gapMask,
                    "overlap_masks": overlapMask,
                    "combined_masks": allMasks,
                    "maximum_overlap": max_overlap,
                    "maximum_gap": max_gap,
                    "channel_percent_available": channel_per,
                }

            elif len(data[k]) == 1:
                try:
                    pct = percent_av[0]
                except TypeError:
                    pct = percent_av
                additions = {
                    "station_completeness": pct,
                    "network": network,
                    "station": station,
                    "location": location,
                    "total_gaps": total_gaps,
                    "total_overlaps": num_overlaps_per_sta,
                    "metric_name": "gapMetricStation",
                    "start_time": trM[0].stats.starttime,
                    "end_time": trM[0].stats.endtime,
                    "gap_masks": gapMask,
                    "overlap_masks": overlapMask,
                    "combined_masks": allMasks,
                    "maximum_overlap": max_overlap,
                    "maximum_gap": max_gap,
                    "channel_percent_available": channel_per,
                }

            else:
                try:
                    if percent_av[k] == 100:
                        additions = {
                            "station_completeness": percent_av[k],
                            "network": network,
                            "station": station,
                            "location": location,
                            "total_gaps": total_gaps,
                            "total_overlaps": total_overlaps,
                            "metric_name": "gapMetricStation",
                            "start_time": trM[0].stats.starttime,
                            "end_time": trM[0].stats.endtime,
                            "gap_masks": gapMask,
                            "overlap_masks": overlapMask,
                            "combined_masks": allMasks,
                            "maximum_overlap": max_overlap,
                            "maximum_gap": max_gap,
                            "channel_percent_available": channel_per,
                        }
                except TypeError:
                    if percent_av == 100:
                        additions = {
                            "station_completeness": percent_av,
                            "network": network,
                            "station": station,
                            "location": location,
                            "total_gaps": total_gaps,
                            "total_overlaps": total_overlaps,
                            "metric_name": "gapMetricStation",
                            "start_time": trM[0].stats.starttime,
                            "end_time": trM[0].stats.endtime,
                            "gap_masks": gapMask,
                            "overlap_masks": overlapMask,
                            "combined_masks": allMasks,
                            "maximum_overlap": max_overlap,
                            "maximum_gap": max_gap,
                            "channel_percent_available": channel_per,
                        }
                else:
                    total_percent_avail = [x for x in percent_av if x != 100]
                    additions = {
                        "station_completeness": np.mean(total_percent_avail),
                        "network": network,
                        "station": station,
                        "location": location,
                        "total_gaps": total_gaps,
                        "total_overlaps": total_overlaps,
                        "metric_name": "gapMetricStation",
                        "start_time": trM[0].stats.starttime,
                        "end_time": trM[0].stats.endtime,
                        "gap_masks": gapMask,
                        "overlap_masks": overlapMask,
                        "combined_masks": allMasks,
                        "maximum_overlap": max_overlap,
                        "maximum_gap": max_gap,
                        "channel_percent_available": channel_per,
                    }
            sum_data.append(additions)

    if database is not None:
        database.insert_metric((sum_data, data))

    return (sum_data, data)
