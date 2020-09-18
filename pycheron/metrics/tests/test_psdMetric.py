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
from pycheron.metrics.psdMetric import psdMetric
import pycheron.metrics.tests.utils as utils


expected_output = {
    "noise2_mask": None,
    "dead_channel_exponent": 0.38710058637617956,
    "percent_above_nhnm": 71.10726643598615,
    "dead_channel_linear": 98.07928228348582,
    "start_time": "2011-01-01T02:16:08.968300",
    "metric_name": "psdMetric",
}


@pytest.mark.parametrize(
    "dead_chan_exp, masks_by_time, by_hour_on",
    [
        (0.5147317990092936, True, True),
        (0.5147317990092936, False, True),
        (0.5147317990092936, True, False),
    ],
)
def test_psdMetric(dead_chan_exp, masks_by_time, by_hour_on):
    stream = utils.get_stream_from_file("test_data/BGU_HHE_002.mseed")
    psd = psdMetric(
        stream,
        byHourOn=by_hour_on,
        generateMasks=True,
        masksByTime=masks_by_time,
    )
    val = [v["dead_channel_exponent"] for v in psd]
    assert utils.percent_diff(dead_chan_exp, val[0], 0.001)
    assert psd[0]["metric_name"] == "psdMetric"


def test_psdMetricEmpty():
    psd = psdMetric([])
    assert psd is None
