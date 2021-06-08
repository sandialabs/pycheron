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
import numpy as np

from pycheron.psd.noise.findOutliers import hampel, hampelFilter, findOutliers


@pytest.mark.parametrize("nwin", [10, 100])
def test_hampel_returns_single_positive_float_value(nwin, file_assets):
    # TODO: Setting nwin to 1 or 2 produces a "invalid value encountered in divide"
    # warning. Unsure if this should be an error. Negative values and zero will
    # raise uncaught exception.
    # ind less than 1/2 the value of nwin will produce uncaught exceptions.
    tr = file_assets["test_traces"][3]
    result = hampel(tr.data, nwin=nwin, ind=50)
    assert isinstance(result, np.float64)
    assert result >= 0


@pytest.mark.parametrize("nwin", [10, 100])
def test_hampelFilter_returns_same_length_array_floats_and_nans(nwin, file_assets):
    # TODO: hampelFilter has a warning that the increment parameter should never
    # be set to anything other than one. If this is true, we should probably not
    # have it as a parameter and just make it a constant or get rid of it entirely.
    # This function also produces "invalid value encountered in greater_equal"
    # for nwin = 10.
    tr = file_assets["test_traces"][3]
    results = hampelFilter(tr.data, nwin=nwin, increment=1)
    assert isinstance(results[0], np.float64)
    assert len(results) == len(tr.data)
    # results at begining and end will be nan becuase of sliding window.
    half_window = nwin // 2
    assert np.isnan(results[: half_window - 1]).all()
    assert np.isnan(results[-half_window:]).all()
    # all other should be positive becuase hampel returns positive nums.
    assert (results[half_window:-half_window] >= 0).all()


def test_findOutliers_default_params_returns_list(file_assets):
    trace = file_assets["test_traces"][3]
    result = findOutliers(trace)
    assert result == []


def test_findOutliers_low_threshold_returns_indicies(file_assets):
    trace = file_assets["test_traces"][3]
    results = findOutliers(trace, threshold=1)
    assert all([result >= 0 for result in results])
    assert all([result < len(trace) for result in results])
    assert all([isinstance(result, int) for result in results])
    assert results[0] == 1131
    assert results[-1] == 206674


def test_findOutliers_low_threshold_with_selectivity_returns_indicies(file_assets):
    trace = file_assets["test_traces"][3]
    results = findOutliers(trace, threshold=1, selectivity=0.1, fixedThreshold=False)
    assert all([result >= 0 for result in results])
    assert all([result < len(trace) for result in results])
    assert all([isinstance(result, int) for result in results])
    assert results[0] == 20
    assert results[-1] == 215979


def test_findOutliers_very_high_threshold_returns_empty_list(file_assets):
    trace = file_assets["test_traces"][3]
    results = findOutliers(trace, threshold=1000000000)
    assert results == []


def test_findOutliers_greater_than_half_of_window_all_same_data_returns_None(
    file_assets,
):
    trace = file_assets["test_traces"][3].copy()
    trace.data = np.full(len(trace.data), fill_value=1)
    result = findOutliers(trace, threshold=1)
    assert result is None
