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

from pycheron.metrics.snrMetric import snrMetric
import pytest
import pycheron.metrics.tests.utils as utils
import obspy.core.event as event
from obspy.core.utcdatetime import UTCDateTime
import numpy as np


"""
When the pick algorithm is used the key for the snr metric is lowercase, snr,
however when the other two algorithms are used the key is uppercase, SNR.
This might cause some issues if the result is passed into other functions

Condition where if (tr.stats.endtime - tr.stats.starttime) < 60:
                    snr = 'Not calculated because data does not fill the window'
However, the algorithm does not exit and continues to evaluate snrMetric. 
Is this behavior what is wanted, or should it exit? 
"""


expected_metric = {
    "SNR": 9.82051194437024,
    "end_time": "2010-02-27T07:16:14.969538",
    "algorithm": "splitWindow",
    "masks": None,
    "windowSecs": 60,
    "start_time": "2010-02-27T06:16:15.019538",
    "snclq": "IU.ANMO.00.BHZ.M",
    "metric_name": "snrMetric",
}

expected_metric_sta = {
    "SNR": 228.66017159000987,
    "end_time": "2011-01-01T23:48:34.158300",
    "algorithm": "staltaTrigger",
    "masks": None,
    "windowSecs": 60,
    "start_time": "2011-01-01T00:00:00.008300",
    "snclq": "UU.ZNPU..HHN.M",
    "metric_name": "snrMetric",
}

expected_metric_pick = {
    "snr": 0.9379343557465231,
    "end_time": "2010-02-27T07:16:14.969538",
    "algorithm": "pick",
    "masks": None,
    "windowSecs": 60,
    "start_time": "2010-02-27T06:16:15.019538",
    "snclq": "IU.ANMO.00.BHZ.M",
    "metric_name": "snrMetric",
}

expected_masks_split_win = {
    "start_time": "2010-02-27T06:16:15.019538",
    "end_time": "2010-02-27T07:16:14.969538",
}
expected_masks_split_stalta = {
    "start_time": "2011-01-01T00:00:00.008300",
    "end_time": "2011-01-01T23:48:34.158300",
}
expected_masks_pick = {
    "start_time": "2010-02-27T06:16:15.019538",
    "end_time": "2010-02-27T07:16:14.969538",
}


@pytest.mark.parametrize(
    "input_type, algorithm_name, expected_output, generate_masks, masks_by_time",
    [
        ("client", "splitWindow", expected_metric, False, True),
        ("file", "staltaTrigger", expected_metric_sta, False, True),
        ("client", "pick", expected_metric_pick, False, True),
        ("client", "splitWindow", expected_masks_split_win, True, True),
        ("file", "staltaTrigger", expected_masks_split_stalta, True, True),
        ("client", "pick", expected_masks_pick, True, True),
        ("client", "splitWindow", expected_masks_split_win, True, False),
        ("file", "staltaTrigger", expected_masks_split_stalta, True, False),
        ("client", "pick", expected_masks_pick, True, False),
    ],
)
def test_snrMetric(input_type, algorithm_name, expected_output, generate_masks, masks_by_time):
    if input_type == "client":
        stream = utils.get_stream_from_client("IU", "ANMO", "00", "BHZ", "2010-02-27T06:16:15.000", 60 * 60)
    else:
        stream = utils.get_stream_from_file("test_data/ZNPU_HHN_001.mseed")
    if algorithm_name == "pick":
        pick = event.origin.Pick(time=UTCDateTime(2010, 2, 27, 6, 18, 8, 968300))
        stream[0].stats.first_arrival = pick
    snr = snrMetric(
        stream,
        algorithm=algorithm_name,
        windowSecs=60,
        generateMasks=generate_masks,
        masksByTime=masks_by_time,
    )
    if generate_masks and masks_by_time:
        utils.compare_dicts(expected_output, snr[0]["masks"][0])
    elif not masks_by_time:
        assert isinstance(snr[0]["masks"], np.ndarray)
    else:
        utils.compare_dicts(expected_output, snr[0])


def test_snrMetricNoAlg():
    stream = utils.get_stream_from_client("IU", "ANMO", "00", "BHZ", "2010-02-27T06:16:15.000", 60 * 60)
    snr = snrMetric(stream, algorithm="test")
    assert snr == []


def test_snrMetricNoData():
    stream = utils.get_stream_from_client("IU", "ANMO", "00", "BHZ", "2010-02-27T06:16:15.000", 30)
    snr = snrMetric(stream)
    assert snr[0]["SNR"] == "Not calculated because data does not fill the window"
