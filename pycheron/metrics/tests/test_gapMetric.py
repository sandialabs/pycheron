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

import pytest
from pycheron.metrics.gapMetric import gapMetric
from obspy.core.utcdatetime import UTCDateTime
import pycheron.metrics.tests.utils as utils
import numpy as np


expected_gap_sum = {
    "overlap_masks": None,
    "station_completeness": 86.73173645157124,
    "total_gaps": 2,
    "start_time": UTCDateTime(2007, 12, 31, 23, 59, 59, 765000),
    "total_overlaps": 3,
    "combined_masks": None,
    "metric_name": "gapMetricStation",
    "network": "BW",
    "maximum_overlap": 2.06000018119812,
    "channel_percent_available": 86.73173645157124,
    "station": "BGLD",
    "end_time": UTCDateTime(2008, 1, 1, 0, 3, 21, 600000),
    "gap_masks": None,
    "maximum_gap": 16.479999780654907,
    "location": "",
}

expected_gap_dets0 = {
    "metric_name": "gapMetricChannel",
    "gap_overlap_start_time": "2008-01-01T00:00:30.620000",
    "gap_overlap_end_time": "2008-01-01T00:00:40.925000",
    "duration": 10.299999952316284,
    "type": "Gap",
    "channel": "EHE",
}

expected_overlap_masks = {
    "end_time": "2007-12-31T23:59:59.765000",
    "start_time": "2007-12-31T23:59:59.760000",
}
expected_gap_masks = {
    "end_time": "2007-12-31T23:59:59.765000",
    "start_time": "2007-12-31T23:59:59.760000",
}


expected_gap_masks_N = {
    "end_time": "2011-01-02T00:00:00.008300",
    "start_time": "2011-01-01T23:59:59.998300",
}
expected_gap_dets_N = {
    "metric_name": "gapMetricChannel",
    "gap_overlap_start_time": "2011-01-02T01:42:43.958300",
    "gap_overlap_end_time": "2011-01-02T01:44:08.978300",
    "duration": 85.00999999046326,
    "type": "Gap",
    "channel": "HHN",
}

expected_gap_masks_Z = {
    "end_time": "2011-01-01T00:00:00.008300",
    "start_time": "2010-12-31T23:59:59.998300",
}
expected_gap_dets_Z = {
    "metric_name": "gapMetricChannel",
    "gap_overlap_start_time": "2011-01-01T02:56:50.958300",
    "gap_overlap_end_time": "2011-01-01T02:59:35.368300",
    "duration": 164.40000009536743,
    "type": "Gap",
    "channel": "HHZ",
}

expected_gap_masks_Z_CTU = {
    "end_time": "2011-01-01T00:00:00.011000",
    "start_time": "2010-12-31T23:59:59.991000",
}
expected_gap_dets_Z_CTU = {
    "metric_name": "gapMetricChannel",
    "gap_overlap_start_time": "2011-01-01T17:48:09.351000",
    "gap_overlap_end_time": "2011-01-01T17:48:09.381000",
    "duration": 0.019999980926513672,
    "type": "Gap",
    "channel": "HHZ",
}

expected_gap_dets_time_jump = {
    "metric_name": "gapMetricChannel",
    "gap_overlap_start_time": "2015-07-01T01:10:20.642500",
    "gap_overlap_end_time": "2015-07-01T01:10:20.027500",
    "duration": -0.6650002002716064,
    "type": "Overlap",
    "channel": "SHZ",
}


def test_gapMetric():
    stream = utils.get_stream_from_file("test_data/sglchan_multgap_multoverlap.mseed")
    gap_sum, gap_dets = gapMetric(stream)
    utils.compare_dicts(expected_gap_sum, gap_sum[0])
    assert len(gap_dets[0]) == 5
    utils.compare_dicts(expected_gap_dets0, gap_dets[0][0])


@pytest.mark.parametrize("masks_by_time, separate_masks", [(True, True), (False, True), (True, False)])
def test_gapMetricMasks(masks_by_time, separate_masks):
    stream = utils.get_stream_from_file("test_data/sglchan_multgap_multoverlap.mseed")
    gap_sum, gap_dets = gapMetric(
        stream,
        generateMasks=True,
        masksByTime=masks_by_time,
        separateMasks=separate_masks,
    )
    if masks_by_time:
        utils.compare_dicts(expected_overlap_masks, gap_sum[0]["overlap_masks"][0])
        utils.compare_dicts(expected_gap_masks, gap_sum[0]["gap_masks"][0])
    else:
        assert isinstance(gap_sum[0]["overlap_masks"], np.ndarray)
        assert isinstance(gap_sum[0]["gap_masks"], np.ndarray)
    utils.compare_dicts(expected_gap_dets0, gap_dets[0][0])


@pytest.mark.parametrize(
    "input_file_name, expected_gap_masks, expected_gap_dets, expected_gap_dets_len, has_gap_masks",
    [
        (
            "test_data/ZNPU_HHN_002.mseed",
            expected_gap_masks_N,
            expected_gap_dets_N,
            7,
            True,
        ),
        (
            "test_data/ZNPU_HHZ_001.mseed",
            expected_gap_masks_Z,
            expected_gap_dets_Z,
            6,
            True,
        ),
        ("test_data/time_jump.mseed", [], expected_gap_dets_time_jump, 1, False),
    ],
)
def test_gapMetricName(
    input_file_name,
    expected_gap_masks,
    expected_gap_dets,
    expected_gap_dets_len,
    has_gap_masks,
):
    stream = utils.get_stream_from_file(input_file_name)
    gap_sum, gap_dets = gapMetric(stream, generateMasks=True)
    if has_gap_masks:
        utils.compare_dicts(expected_gap_masks, gap_sum[0]["gap_masks"][0])
    else:
        assert gap_sum[0]["gap_masks"] == []
    assert len(gap_dets[0]) == expected_gap_dets_len
    utils.compare_dicts(expected_gap_dets, gap_dets[0][0])
