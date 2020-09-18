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

from pycheron.metrics.repeatedAmplitude import repeatedAmplitudeMetric
import pycheron.metrics.tests.utils as utils

expected_output = {
    "count": 0,
    "metric_name": "repeatedAmplitudeMetric",
    "rep_amp": [],
    "start_time": "2010-02-27T06:16:15.019538",
    "mask": None,
    "snclq": "IU.ANMO.00.BHZ",
    "end_time": "2010-02-28T06:16:14.969538",
}
expected_rep_amp = {
    "duration": 0.2,
    "start_time": "2012-12-12T07:45:29.468400",
    "end_time": "2012-12-12T07:45:29.648400",
    "value": "1292",
}


def test_repeatedAmplitudes():
    st0 = utils.get_stream_from_client(
        "AK", "GHO", "", "BHN", "2012-12-12T00:00:00.000", 60 * 60 * 24
    )
    st = utils.get_stream_from_client(
        "IU", "ANMO", "00", "BHZ", "2010-02-27T06:16:15.000", 60 * 60 * 24
    )
    repAmp = repeatedAmplitudeMetric(st, 10)
    assert len(repAmp) == 1
    utils.compare_dicts(expected_output, repAmp[0])

    repAmp = repeatedAmplitudeMetric(st0, 10)
    assert repAmp[0]["count"] == 1
    utils.compare_dicts(expected_rep_amp, repAmp[0]["rep_amp"][0])

    repAmp = repeatedAmplitudeMetric(st0, 10, generateMasks=True)
    assert repAmp[0]["count"] == 1
    utils.compare_dicts(expected_rep_amp, repAmp[0]["rep_amp"][0])

    st0 = utils.get_stream_from_client(
        "AK", "GHO", "", "BHN", "2012-12-12T00:00:00.000", 60 * 60 * 24
    )
    repAmp = repeatedAmplitudeMetric(st0, 5, generateMasks=True)
    assert repAmp[0]["count"] == 180
