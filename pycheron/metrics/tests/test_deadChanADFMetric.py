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
from pycheron.metrics.deadChanADFMetric import deadChanADFMetric
import pycheron.metrics.tests.utils as utils


expected_metric1_mbt_true = {
    "metric_name": "deadChanADFMetric",
    "snclq": "BW.BGLD..EHE",
    "start_time": "2007-12-31T23:59:59.765000",
    "masks": "No Masks Created",
    "end_time": "2008-01-01T00:00:30.620000",
    "dead_channel_times": None,
}

expected_metric1_mbt_false = {
    "metric_name": "deadChanADFMetric",
    "snclq": "BW.BGLD..EHE",
    "start_time": "2007-12-31T23:59:59.765000",
    "masks": "No Masks Created",
    "end_time": "2008-01-01T00:00:30.620000",
    "dead_channel_times": None,
}

expected_metric1_sampl_rate100 = {
    "metric_name": "deadChanADFMetric",
    "snclq": "BW.BGLD..EHE",
    "start_time": "2007-12-31T23:59:59.765000",
    "masks": "No Masks Created",
    "end_time": "2008-01-01T00:01:01.475000",
    "dead_channel_times": None,
}

expected_metric1_sampl_rate40 = {
    "metric_name": "deadChanADFMetric",
    "snclq": "BW.BGLD..EHE",
    "start_time": "2007-12-31T23:59:59.765000",
    "masks": "No Masks Created",
    "end_time": "2008-01-01T00:02:34.040000",
    "dead_channel_times": None,
}

expected_metric1_sampl_rate20 = {
    "metric_name": "deadChanADFMetric",
    "snclq": "BW.BGLD..EHE",
    "start_time": "2007-12-31T23:59:59.765000",
    "masks": "No Masks Created",
    "end_time": "2008-01-01T00:05:08.315000",
    "dead_channel_times": None,
}


@pytest.mark.parametrize(
    "expected_output, masks_by_time, win_size, gen_masks",
    [
        (expected_metric1_mbt_true, True, 500, True),
        (expected_metric1_mbt_false, False, 500, True),
        (expected_metric1_mbt_true, True, None, True),
        (expected_metric1_mbt_true, True, None, False),
    ],
)
def test_deadChanADFMetric(expected_output, masks_by_time, win_size, gen_masks):
    stream = utils.get_stream_from_file("test_data/sglchan_multgap_multoverlap.mseed")
    metric = deadChanADFMetric(
        stream, generateMasks=gen_masks, masksByTime=masks_by_time, win_size=win_size
    )
    utils.compare_dicts(expected_output, metric[0])


@pytest.mark.parametrize(
    "expected_output, sample_rate, p_value, masks_by_time",
    [
        (expected_metric1_sampl_rate100, 100, "p_value", True),
        (expected_metric1_sampl_rate40, 40, "p_value", True),
        (expected_metric1_sampl_rate20, 20, "p_value", True),
        (expected_metric1_sampl_rate100, 100, None, False),
    ],
)
def test_deadChanADFMetricSampleRate(
    expected_output, sample_rate, p_value, masks_by_time
):
    stream = utils.get_stream_from_file("test_data/sglchan_multgap_multoverlap.mseed")
    for i in range(len(stream)):
        stream[i].stats.sampling_rate = sample_rate
    metric = deadChanADFMetric(
        stream,
        generateMasks=True,
        win_size=None,
        use_thresh=p_value,
        masksByTime=masks_by_time,
    )
    utils.compare_dicts(expected_output, metric[0])
