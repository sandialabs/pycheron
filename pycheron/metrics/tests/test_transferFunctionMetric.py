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

from pycheron.metrics.transferFunctionMetric import transferFunctionMetric
import pycheron.metrics.tests.utils as utils


expected_output = {
    "metric_name": "transferFunctionMetric",
    "start_time": "2011-05-01T01:00:00.019500",
    "snclq": "CI.PASC.00.BHZ:CI.PASC.10.BHZ",
    "end_time": "2011-05-01T01:59:59.994500",
    "ms_coherence": 0.9999875790954714,
    "gain_ratio": 0.9058981533571073,
    "phase_diff": 0.13553317118827302,
}

expected_output_update = {
    "metric_name": "transferFunctionMetric",
    "start_time": "2011-05-01T01:00:00.019500",
    "snclq": "CI.PASC.00.BHZ:CI.PASC.10.BHZ",
    "end_time": "2011-05-01T01:59:59.994500",
    "ms_coherence": 0.37147995324408395,
    "gain_ratio": 0.29845579397758604,
    "phase_diff": -0.300582352310399,
}


def test_transferFunctionMetric():
    st1 = utils.get_stream_from_client("CI", "PASC", "00", "BHZ", "2011-05-01T01:00:00.000", 3600)
    st2 = utils.get_stream_from_client("CI", "PASC", "10", "BHZ", "2011-05-01T01:00:00.000", 3600)
    tr1 = st1[0]
    tr2 = st2[0]
    tf_metric = transferFunctionMetric(tr1, tr2)
    utils.compare_dicts(expected_output, tf_metric)
