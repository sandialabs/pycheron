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
from pycheron.rollseis.roll_stalta import roll_stalta, stalta_python


@pytest.fixture
def stalta_assets():
    return {"x": np.repeat([1, 5, 3, 2, 1], 20), "n_sta": 3, "n_lta": 6, "increment": 1}


def test_n_sta_gt_len_roll_stalta(stalta_assets):
    n_sta = len(stalta_assets["x"]) + 1
    stalta = roll_stalta(
        stalta_assets["x"], n_sta, stalta_assets["n_lta"], stalta_assets["increment"]
    )
    assert stalta is None


def test_n_lta_gt_len_roll_stalta(stalta_assets):
    n_lta = len(stalta_assets["x"]) + 1
    stalta = roll_stalta(
        stalta_assets["x"], stalta_assets["n_sta"], n_lta, stalta_assets["increment"]
    )
    assert stalta is None


@pytest.mark.parametrize("increment", list(range(-10, 1)))
def test_increment_lt_1_roll_stalta(stalta_assets, increment):
    stalta = roll_stalta(
        stalta_assets["x"], stalta_assets["n_sta"], stalta_assets["n_lta"], increment
    )
    assert stalta is None


def test_valid_roll_stalta(stalta_assets):
    stalta = roll_stalta(
        stalta_assets["x"],
        stalta_assets["n_sta"],
        stalta_assets["n_lta"],
        stalta_assets["increment"],
    )
    grab_expected_reals = [
        isinstance(x, np.float64)
        for x in stalta[stalta_assets["n_lta"] : -stalta_assets["n_sta"]]
    ]
    grab_expected_nans_lta = [np.isnan(x) for x in stalta[: stalta_assets["n_lta"]]]
    grab_expected_nans_sta = [np.isnan(x) for x in stalta[-stalta_assets["n_sta"] :]]

    assert isinstance(stalta, np.ndarray)
    assert len(stalta) > 0
    assert len(stalta) == len(stalta_assets["x"])
    assert all(x for x in grab_expected_reals)
    assert all(not np.isnan(x) for x in grab_expected_reals)
    assert all(x for x in grab_expected_nans_lta)
    assert all(x for x in grab_expected_nans_sta)


@pytest.mark.parametrize("n_sta", list(range(1, 100)))
def test_num_nans_n_sta_roll_stalta(stalta_assets, n_sta):
    stalta = roll_stalta(
        stalta_assets["x"], n_sta, stalta_assets["n_lta"], stalta_assets["increment"]
    )
    grab_nans = [np.isnan(x) for x in stalta[-n_sta:]]
    assert all(x for x in grab_nans)


def test_n_sta_equals_0_roll_stalta(stalta_assets):
    # TODO: in refactor, ensure roll_stalta has a check for n_sta = 0
    n_sta = 0
    with pytest.raises(ZeroDivisionError):
        stalta = roll_stalta(
            stalta_assets["x"],
            n_sta,
            stalta_assets["n_lta"],
            stalta_assets["increment"],
        )


@pytest.mark.parametrize("n_lta", list(range(1, 100)))
def test_num_nans_n_lta_roll_stalta(stalta_assets, n_lta):
    stalta = roll_stalta(
        stalta_assets["x"], stalta_assets["n_sta"], n_lta, stalta_assets["increment"]
    )
    grab_nans = [np.isnan(x) for x in stalta[:n_lta]]
    assert all(x for x in grab_nans)


def test_n_Lta_equals_0_roll_stalta(stalta_assets):
    # TODO: in refactor, ensure roll_stalta has a check for n_sta = 0
    n_lta = 0
    with pytest.raises(ZeroDivisionError):
        stalta = roll_stalta(
            stalta_assets["x"],
            stalta_assets["n_sta"],
            n_lta,
            stalta_assets["increment"],
        )
