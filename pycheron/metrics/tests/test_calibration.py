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

from pycheron.metrics.calibration import calibrationMetric
import pycheron.metrics.tests.utils as utils
import pytest

expected_detection0 = {
    "env_median": 104.05356500003606,
    "cal_start_time": "2011-01-06T02:10:54.388300",
    "cal_end_time": "2011-01-06T02:11:32.868300",
    "max_ySmooth": 6.081637831289797,
    "max_stalta": 604.1447908022914,
    "stft_var": 2.181164305836342,
    "env_mean": 174.0817802294856,
    "env_var": 15399.731750251214,
    "initial_cal_start_time": "2011-01-06T02:10:45.878300",
    "initial_cal_end_time": "2011-01-06T02:12:19.758300",
}

expected_output_empty = {
    "metric_name": "calibrationMetric",
    "start_time": "2011-01-01T02:16:08.968300",
    "snclq": "UU.BGU..HHE",
    "detections": [],
    "num_cals_detected": 0,
    "end_time": "2011-01-01T20:11:23.328300",
}


@pytest.mark.parametrize(
    "input_file_name, expected_output",
    [
        ("test_data/ZNPU_HHN_006.mseed", expected_detection0),
        ("test_data/BGU_HHE_001.mseed", expected_output_empty),
    ],
)
def test_calibration(input_file_name, expected_output):
    stream = utils.get_stream_from_file(input_file_name)
    calibration = calibrationMetric(stream)
    if calibration[0]["detections"]:
        assert len(calibration[0]["detections"]) == 2
        utils.compare_dicts(expected_output, calibration[0]["detections"][0])
    else:
        print(expected_output)
        print(calibration[0])
        utils.compare_dicts(expected_output, calibration[0])
