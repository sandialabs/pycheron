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

import numpy as np
import pytest

from pycheron.sigpro.STALTA import STALTA


@pytest.fixture(scope="module")
def short_trace(file_assets):
    test_trace_copy = file_assets["test_trace"].copy()
    test_trace_copy.data = test_trace_copy.data[:50000]
    test_trace_copy.stats["npts"] = 50000
    return test_trace_copy


@pytest.fixture(scope="module")
def short_default_STALTA(short_trace):
    stalta_applied_trace = STALTA(short_trace)
    return stalta_applied_trace


def assert_is_np_array_of_floats(array_to_test):
    assert isinstance(array_to_test, np.ndarray)
    assert isinstance(array_to_test[0], np.float64)


def test_STALTA_returns_same_length_array(short_trace, short_default_STALTA):
    # trim trace so it is quicker to run.
    original_trace_length = len(short_trace)
    new_STALTA_length = len(short_default_STALTA)
    assert original_trace_length == new_STALTA_length


def test_STALTA_default_params_returns_float_np_array(short_default_STALTA):
    assert_is_np_array_of_floats(short_default_STALTA)


def test_STALTA_returns_nans_for_edges_of_trace(short_default_STALTA):
    # not enough data available for sliding window algorithms in data near
    # begining or end of trace. Nans are returned in these positions.
    assert all(np.isnan(x) for x in short_default_STALTA[:600])
    assert all(np.isnan(x) for x in short_default_STALTA[-50:])


def test_increasing_increment_param_returns_more_nan_values(
    short_trace, short_default_STALTA
):
    # increment determines how many values to skip between testing points.
    # default is 1 so it tests every data point.
    default_STALTA_num_nans = np.sum(np.isnan(short_default_STALTA))

    increment_increased_STALTA = STALTA(short_trace, increment=2)
    increment_increased_num_nans = np.sum(np.isnan(increment_increased_STALTA))

    assert increment_increased_num_nans > default_STALTA_num_nans


def test_setting_demean_and_detrend_params_changes_outputs(
    short_trace, short_default_STALTA
):
    both_demean_and_detrend = short_default_STALTA
    demean_only = STALTA(short_trace, detrend=False)
    detrend_only = STALTA(short_trace, demean=False)
    neither_detrend_nor_demean = STALTA(short_trace, demean=False, detrend=False)

    assert (both_demean_and_detrend != demean_only).any()
    assert (both_demean_and_detrend != detrend_only).any()
    assert (both_demean_and_detrend != neither_detrend_nor_demean).any()

    assert (demean_only != detrend_only).any()
    assert (demean_only != neither_detrend_nor_demean).any()

    assert (detrend_only != neither_detrend_nor_demean).any()


def test_differant_algorithms_produce_differant_results(
    short_trace, short_default_STALTA
):
    classic_LR_algorithm = short_default_STALTA

    classic_RR_algorithm = STALTA(short_trace, algorithm="classic_RR")
    assert_is_np_array_of_floats(classic_RR_algorithm)

    e_and_s_envelope_algorithm = STALTA(
        short_trace, algorithm="EarleAndShearer_envelope"
    )
    assert_is_np_array_of_floats(e_and_s_envelope_algorithm)

    assert (classic_LR_algorithm != classic_RR_algorithm).any()
    assert (classic_LR_algorithm != e_and_s_envelope_algorithm).any()
    assert (classic_RR_algorithm != e_and_s_envelope_algorithm).any()


def test_setting_sta_secs_too_short_returns_none(short_trace):
    # Currently if the staSecs parameter is too small, it logs an error, then
    # returns None. This should probably be changed in the future to raise a
    # ValueError.
    assert STALTA(short_trace, staSecs=0) is None


def test_too_short_trace_returns_none(short_trace):
    # Currently if the trace is shorter than the ltaSecs parameter
    # multiplied by the sampling rate, it logs an error and
    # returns None. This should probably be changed in the future to raise a
    # ValueError.
    invalidly_short_trace = short_trace.copy()
    invalidly_short_trace.data = invalidly_short_trace.data[:100]
    invalidly_short_trace.stats["npts"] = 100

    assert STALTA(invalidly_short_trace, ltaSecs=100) is None


def test_unknown_algorithm_returns_none(short_trace):
    # Currently if the algorithm parameter is not recognized, it logs an error
    # and returns None. This should probably be changed in the future to raise a
    # ValueError.
    assert STALTA(short_trace, algorithm="dijkstras") is None


def test_E_and_S_algorithm_with_increment_greater_than_1_returns_none(short_trace):
    # Currently if the algorithm parameter is not recognized, it logs an error
    # and returns None. This should probably be changed in the future to raise a
    # ValueError.
    result = STALTA(short_trace, algorithm="EarleAndShearer_envelope", increment=2)
    assert result is None


def test_RR_algorithm_with_increment_greater_than_1_returns_none(short_trace):
    # Currently if the algorithm parameter is not recognized, it logs an error
    # and returns None. This should probably be changed in the future to raise a
    # ValueError.
    result = STALTA(short_trace, algorithm="classic_RR", increment=2)
    assert result is None
