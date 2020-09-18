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

from pycheron.metrics.basicStatsMetric import basicStatsMetric
from obspy.core.utcdatetime import UTCDateTime
import pytest
import pycheron.metrics.tests.utils as utils


expected_output = {
    "rms": 866.4443879572756,
    "number_of_samples": 6451437,
    "metric_name": "basicStatsMetric",
    "max_mask": None,
    "std_mask": None,
    "standard_deviation": 197.98584010776364,
    "start_time": UTCDateTime(2011, 1, 1, 2, 16, 8, 968300),
    "median": 855.0,
    "maximum": 9223,
    "snclq": "UU.BGU..HHE",
    "min_mask": None,
    "variance_mask": None,
    "minimum": -9569,
    "end_time": UTCDateTime(2011, 1, 1, 20, 11, 23, 328300),
    "median_mask": None,
    "variance": 39198.39288317695,
    "mean_mask": None,
    "rms_mask": None,
    "mean": 843.5208856569475,
}


@pytest.mark.parametrize(
    "masks_by_time, threshold", [(True, None), (False, None), (False, 10)]
)
def test_basicStatsMetric(masks_by_time, threshold):
    stream = utils.get_stream_from_file("test_data/BGU_HHE_001.mseed")
    basic_stats = basicStatsMetric(
        stream, generateMasks=True, masksByTime=masks_by_time, maxThreshold=threshold
    )
    if threshold:
        assert len(basic_stats[0]["max_mask"]) > 0
    else:
        utils.compare_dicts(expected_output, basic_stats[0])
