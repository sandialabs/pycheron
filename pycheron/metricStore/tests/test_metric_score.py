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

from pycheron.metricStore.metric_store import MetricStore, MultiSnclqMetricStore
from obspy import UTCDateTime
import numpy as np
from pycheron.metrics.tests.utils import percent_diff


expected_vals = [-0.03413873, -0.04130493, -0.20611203]


class YearMonth(object):
    def __init__(self, year, month):
        self.year = year
        self.month = month


def add_data_to_metricStore(ms):
    snclq = "LPAZ.BHZ"
    metric_name = "calibration.env_var"
    ms.metric_names(snclq)

    for number_of_days in range(60):
        # move a day forward
        date = UTCDateTime(2012, 8, 28) + (24 * 3600 * number_of_days)
        metric_value = np.random.randn(3)
        ms.put(snclq, metric_name, metric_value, date, overwrite=True, reshape=False)

    return ms.get_days(snclq, metric_name)


def test_value_props_single():
    metric_store_single = MetricStore(id=0)
    values, dates = add_data_to_metricStore(metric_store_single)
    assert isinstance(values, list)
    assert len(values) > 0


def test_value_props_multi():
    metric_store_mult = MultiSnclqMetricStore(id=0)
    values, dates = add_data_to_metricStore(metric_store_mult)
    assert isinstance(values, list)
    assert len(values) > 0


def test_year_month_index():
    metric_store_mult = MultiSnclqMetricStore(id=0)
    metric_store_single = MetricStore(id=0)

    ym = YearMonth(2012, 12)
    ym_single = metric_store_single._year_month_index(ym)
    ym_mult = metric_store_mult._year_month_index(ym)
    assert ym_single == "2012-12"
    assert ym_mult == "2012-12"


def test_metric_store_single_key_error():
    metric_store_single = MetricStore(id=0)
    snclq = "LPAZ.BHZ"
    metric_name = "fake_metric.env_var"
    add_data_to_metricStore(metric_store_single)
    data = metric_store_single.get(snclq, metric_name, (2012, 9))
    assert data is None


def test_metric_store_multi_key_error():
    metric_store_mult = MultiSnclqMetricStore(id=0)
    snclq = "LPAZ.BHZ"
    metric_name = "fake_metric.env_var"
    add_data_to_metricStore(metric_store_mult)
    data = metric_store_mult.get(snclq, metric_name, (2012, 9))
    assert data is None


def test_metric_store_multi():
    metric_store_mult = MultiSnclqMetricStore(id=0)
    snclq = "LPAZ.BHZ"
    metric_name = "calibration.env_var"
    values, dates = add_data_to_metricStore(metric_store_mult)
    assert len(values) == 31
    assert len(values[0]) == 3

    val = metric_store_mult.get(snclq, metric_name, (2012, 9))
    for i in range(3):
        assert percent_diff(expected_vals[i], val[i], 0.0001)
    data = metric_store_mult._day_tracker(snclq)
    assert isinstance(data, MultiSnclqMetricStore._DayTracker)


def test_metric_store_single():
    metric_store_single = MetricStore(id=0)
    snclq = "LPAZ.BHZ"
    metric_name = "calibration.env_var"
    values, dates = add_data_to_metricStore(metric_store_single)
    assert len(values) == 31
    assert len(values[0]) == 3

    val = metric_store_single.get(snclq, metric_name, (2012, 9))
    for i in range(3):
        assert percent_diff(expected_vals[i], val[i], 0.0001)
    data = metric_store_single._day_tracker.get_year_month_data()
    assert data == (None, None)
