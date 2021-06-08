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

from pycheron.metrics.correlationMetric import correlationMetric
import pytest
from obspy.core.utcdatetime import UTCDateTime
import pycheron.metrics.tests.utils as utils

expected_value = {
    "p_value": 0.0,
    "metric_name": "correlationMetric",
    "start_time": UTCDateTime(2013, 3, 1, 0, 0, 0, 19500),
    "correlation_coefficient": -0.2041723553359387,
    "snclq": "IU.ANMO.00.BH1.M:IU.ANMO.00.BH2.M",
    "end_time": UTCDateTime(2013, 3, 1, 23, 59, 59, 969500),
}


@pytest.mark.parametrize("input_file_name", ["test_data/ZNPU_HHN_001.mseed", "test_data/CTU_HHE_004.mseed"])
def test_correlationMetricNone(input_file_name):
    stream = utils.get_stream_from_file(input_file_name)
    corr = []
    for i in range(len(stream) - 1):
        corr.append(correlationMetric(stream[i], stream[i + 1]))
    assert not (all(corr))


def test_correlationMetrticValue():
    st = utils.get_stream_from_client("IU", "ANMO", "00", "BH*", "2013-03-01T00:00:00.000", 1440 * 60)
    corr = correlationMetric(st[0], st[1])
    utils.compare_dicts(expected_value, corr)


def test_correlationMetricUnequalEndtime():
    st1 = utils.get_stream_from_client("IU", "ANMO", "00", "BH*", "2013-03-01T00:00:00.000", 1440 * 60)
    st2 = utils.get_stream_from_client("IU", "ANMO", "00", "BH*", "2013-03-01T00:00:00.000", 60 * 60)
    corr = correlationMetric(st1[0], st2[0])
    assert corr is None


def test_correlationMetricLength():
    st = utils.get_stream_from_client("IU", "ANMO", "00", "BH*", "2013-03-01T00:00:00.000", 1440 * 60)
    tr = st[0]
    tr1 = st[1]
    tr.data = tr.data[:-100]
    corr = correlationMetric(tr, tr1)
    assert corr is None


def test_correlationMetricSampling():
    st = utils.get_stream_from_client("IU", "ANMO", "00", "BH*", "2013-03-01T00:00:00.000", 1440 * 60)
    tr = st[0]
    tr1 = st[1]
    tr.stats.sampling_rate = 20
    tr1.stats.sampling_rate = 100
    corr = correlationMetric(tr, tr1)
    assert corr is None


def test_correlationMetricDifStation():
    st = utils.get_stream_from_client("IU", "ANMO", "00", "BH*", "2013-03-01T00:00:00.000", 1440 * 60)
    tr = st[0]
    tr1 = st[1]
    tr.stats.station = "test"
    corr = correlationMetric(tr, tr1)
    assert corr is None
