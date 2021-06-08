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

from pycheron.metrics.staltaMetric import staltaMetric
import pycheron.metrics.tests.utils as utils

expected_output_single = {
    "metric_name": "staltaMetric",
    "event_time": "2012-12-12T21:37:26.828400",
    "start_time": "2012-12-12T00:00:00.008400",
    "snclq": "AK.GHO..BHN.M",
    "max_stalta": 19846,
    "end_time": "2012-12-12T23:59:59.988400",
}

expected_output_multiple = {
    "metric_name": "staltaMetric",
    "event_time": "2008-01-01T00:01:34.485000",
    "start_time": "2007-12-31T23:59:59.765000",
    "snclq": "BW.BGLD..EHE.D",
    "max_stalta": 10720.887,
    "end_time": "2008-01-01T00:03:21.600000",
}


def test_staltaMetricNonPool():
    st = utils.get_stream_from_client("AK", "GHO", "", "BHN", "2012-12-12T00:00:00.000", 60 * 60 * 24)
    stalta = staltaMetric(st)
    utils.compare_dicts(expected_output_single, stalta[0])


def test_staltaMetricPool():
    st1 = utils.get_stream_from_file("test_data/sglchan_multgap_multoverlap.mseed")
    stalta = staltaMetric(st1, processes=2)
    utils.compare_dicts(expected_output_multiple, stalta[0])
