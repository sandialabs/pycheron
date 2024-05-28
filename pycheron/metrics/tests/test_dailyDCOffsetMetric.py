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

from pycheron.metrics.dailyDCOffsetMetric import dailyDCOffSetMetric
from pycheron.metrics.basicStatsMetric import basicStatsMetric
import pycheron.metrics.tests.utils as utils
from obspy.core.utcdatetime import UTCDateTime


expected_output = {
    "daily_dc_offset_value": 0.0025,
    "start_time": UTCDateTime(2007, 12, 31, 23, 59, 59, 915000),
    "end_time": UTCDateTime(2008, 1, 1, 0, 0, 1, 970000),
    "snclq": "BW.BGLD..EHE",
    "metric_name": "dailyDCOffset",
}


def test_dailyDCOffsetMetric():
    stream = utils.get_stream_from_file("test_data/BGU_HHE_001.mseed")
    basic_stats = basicStatsMetric(stream, generateMasks=True)
    offset_metric = dailyDCOffSetMetric(basic_stats, offsetDays=1)
    assert offset_metric[0]["daily_dc_offset_value"] is None


def test_dailyDCOffsetMetricEmpty():
    offset_metric = dailyDCOffSetMetric({})
    assert offset_metric is None


def test_dailyDCOffsetMetricNoMean():
    offset_metric = dailyDCOffSetMetric({"val": [12334, 1233]})
    assert offset_metric is None


def test_DC():
    d = utils.make_dict_of_lists(basicStatsMetric(utils.get_stream_from_file("test_data/qualityflags.mseed")))
    print(d)
    off = dailyDCOffSetMetric(d, offsetDays=1)
    assert len(off) > 0
    utils.compare_dicts(expected_output, off)
    off = dailyDCOffSetMetric(d, offsetDays=1, OutputType=0)
    utils.compare_dicts(expected_output, off[0])
