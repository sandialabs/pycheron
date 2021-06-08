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
from pycheron.metrics.deadChanMeanMetric import deadChanMean
import pycheron.metrics.tests.utils as utils
import numpy as np

masks1 = {
    "Startime": "2008-01-01T00:00:40.925000",
    "Endtime": "2008-01-01T00:01:08.925000",
}

dead_channel_times1 = {
    "dead_starttime": "2008-01-01T00:00:40.925000",
    "dead_endtime": "2008-01-01T00:01:08.925000",
}

metric_output100 = {
    "metric_name": "deadChanMeanMetric",
    "snclq": "BW.BGLD..EHE",
    "start_time": "2007-12-31T23:59:59.765000",
    "masks": None,
    "end_time": "2008-01-01T00:01:01.475000",
    "dead_channel_times": None,
}

metric_output40 = {
    "metric_name": "deadChanMeanMetric",
    "snclq": "BW.BGLD..EHE",
    "start_time": "2007-12-31T23:59:59.765000",
    "masks": None,
    "end_time": "2008-01-01T00:02:34.040000",
    "dead_channel_times": None,
}

metric_output20 = {
    "metric_name": "deadChanMeanMetric",
    "snclq": "BW.BGLD..EHE",
    "start_time": "2007-12-31T23:59:59.765000",
    "masks": None,
    "end_time": "2008-01-01T00:05:08.315000",
    "dead_channel_times": None,
}


@pytest.mark.parametrize(
    "sample_rate, win_size, masks_by_time, is_mask_array",
    [(100.0, 2800, True, False), (100.0, 2800, False, True)],
)
def test_deadChannelMeanMetric(sample_rate, win_size, masks_by_time, is_mask_array):
    stream = utils.get_stream_from_file("test_data/sglchan_multgap_multoverlap.mseed")
    for i in range(len(stream)):
        stream[i].stats.sampling_rate = sample_rate
    metric = deadChanMean(
        stream,
        generateMasks=True,
        threshold=10.0,
        win_size=win_size,
        masksByTime=masks_by_time,
    )
    if is_mask_array:
        assert isinstance(metric[1]["masks"], np.ndarray)
    else:
        utils.compare_dicts(masks1, metric[1]["masks"][0])
    utils.compare_dicts(dead_channel_times1, metric[1]["dead_channel_times"][0])


@pytest.mark.parametrize(
    "expected_output, sample_rate, win_size, masks_by_time",
    [
        (metric_output100, 100.0, None, True),
        (metric_output40, 40.0, None, True),
        (metric_output20, 20.0, None, True),
    ],
)
def test_deadChannelMeanMetricSampRate(expected_output, sample_rate, win_size, masks_by_time):
    stream = utils.get_stream_from_file("test_data/sglchan_multgap_multoverlap.mseed")
    for i in range(len(stream)):
        stream[i].stats.sampling_rate = sample_rate
    metric = deadChanMean(
        stream,
        generateMasks=True,
        threshold=10.0,
        win_size=win_size,
        masksByTime=masks_by_time,
    )
    utils.compare_dicts(expected_output, metric[0])
