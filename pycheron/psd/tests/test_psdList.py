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
from obspy.core.utcdatetime import UTCDateTime
import pytest

from pycheron.psd.psdList import psdList


def test_psdList_returns_None_if_stream_is_less_than_an_hour(file_assets):
    # Unclear if this is an error that should raise a ValueError, or
    # if we really want it returning none.
    stream_less_than_one_hour = file_assets["test_streams"][0]
    results = psdList(stream_less_than_one_hour)
    assert results[0] is None


def test_psdList_returns_correct_num_traces_and_segments(file_assets):
    streams = file_assets["test_streams"][3:7]
    for stream in streams:
        psdList_results = psdList(stream)
        for index, trace in enumerate(stream):
            channel_first_char = trace.stats.channel[0]
            z = 3600
            if channel_first_char == "V":
                z = 24 * 3600
            elif channel_first_char == "L":
                z = 3 * 3600
            elif channel_first_char == "M":
                z = 2 * 3600
            num_segments = 0
            secs_left = trace.stats.endtime - trace.stats.starttime
            while secs_left >= 0.99 * z:
                num_segments += 1
                secs_left -= z / 2

            assert len(psdList_results[index]) == num_segments


def test_psdList_returns_correct_metadata(BHZ_channel_stream_and_psd_list):
    # Bhz channel trace has correct snlcq (station network location channel quality)
    # string and segments have correct start and end times after psdList.
    stream = BHZ_channel_stream_and_psd_list[0]
    psdList_results = BHZ_channel_stream_and_psd_list[1]

    for index, tr in enumerate(stream):
        snlcq_tr = [
            "{0}.{1}.{2}.{3}.{4}".format(
                str(tr.stats.network),
                str(tr.stats.station),
                str(tr.stats.location),
                str(tr.stats.channel),
                str(tr.stats.mseed.dataquality),
            )
        ]
        start_tr = tr.stats["starttime"]
        end_tr = tr.stats["endtime"]
        for index_seg, segment in enumerate(psdList_results[index]):
            assert segment[2] == snlcq_tr
            start_expected = start_tr + 1800 * index_seg
            end_expected = start_expected + 3600
            assert UTCDateTime(segment[3]) == start_expected
            if end_expected > end_tr:
                end_expected = end_tr
            assert UTCDateTime(segment[4]) == end_expected


def test_psdList_returns_correct_types(BHZ_channel_stream_and_psd_list):
    # Testing return types for everything that isn't covered in metadata.
    stream = BHZ_channel_stream_and_psd_list[0]
    psdList_results = BHZ_channel_stream_and_psd_list[1]

    for index, tr in enumerate(stream):
        sum_trace = 0
        for segment in psdList_results[index]:
            # Noise and Frequency Matricies should be floating point numpy arrays
            # of same length.
            assert isinstance(segment[0], np.ndarray)
            assert isinstance(segment[0][0], np.float64)
            assert isinstance(segment[1], np.ndarray)
            assert isinstance(segment[1][0], np.float64)
            assert len(segment[0]) == len(segment[1])
