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
from math import floor, ceil
from pycheron.rollseis.roll_mean import roll_mean


@pytest.mark.parametrize("x", [np.array([1, 2, 3]), []])
@pytest.mark.parametrize("align", ["left", "right", "center"])
def test_nwin_gt_x_roll_mean(roll_assets, x, align):
    mean_vals = roll_mean(x, roll_assets["nwin"], roll_assets["ind"], align)
    assert mean_vals is None


@pytest.mark.parametrize("ind", list(range(-10, 1)))
@pytest.mark.parametrize("align", ["left", "right", "center"])
def test_ind_lt_1_roll_mean(roll_assets, ind, align):
    mean_vals = roll_mean(roll_assets["x"], roll_assets["nwin"], ind, align)
    assert mean_vals is None


@pytest.mark.parametrize("nwin", [6, 7])
@pytest.mark.parametrize("ind", [1, 2, 3, 4, 5])
def test_nans_mean_roll(roll_assets, nwin, ind):
    # TODO: for refactor, ensure "left", "right", and "center" are the only values accepted
    # currently, "right" will produce the same result as a random string
    #
    # For refactor:
    # There needs to be a check for when nwin==0 in the original code
    # This hasn't been implemented but will fail if this is the case
    left = roll_mean(roll_assets["x"], nwin, ind, "left")
    right = roll_mean(roll_assets["x"], nwin, ind, "right")
    center = roll_mean(roll_assets["x"], nwin, ind, "center")

    # because of the rolling window,
    # for a left window: the last nwin-1 values should be NaN
    # for a right window: the first nwin-1 values should be NaN
    # for a center window: the first nwin-1 / 2 and last nwin-1 / 2 values should be NaN
    # (nwin increased by 1 in order to guarantee odd window size if necessary)
    left_nan_window = nwin * -1 + 1
    right_nan_window = nwin - 1

    # determined by size of rolling window
    left_grab_nans = [np.isnan(x) for x in left[left_nan_window:]]
    right_grab_nans = [np.isnan(x) for x in right[:right_nan_window]]

    assert len(left_grab_nans) == len(right_grab_nans)

    center_grab_nans_l = [np.isnan(x) for x in center[int(floor(left_nan_window / 2.0)) :]]
    center_grab_nans_r = [np.isnan(x) for x in center[: int(ceil(right_nan_window / 2.0))]]

    assert len(center_grab_nans_l) == len(center_grab_nans_r)

    assert all(x for x in left_grab_nans)
    assert all(x for x in right_grab_nans)
    assert all(x for x in center_grab_nans_l)
    assert all(x for x in center_grab_nans_r)


@pytest.mark.parametrize("nwin", [6, 7])
def test_reals_mean_roll(roll_assets, nwin):
    left = roll_mean(roll_assets["x"], nwin, roll_assets["ind"], "left")
    right = roll_mean(roll_assets["x"], nwin, roll_assets["ind"], "right")
    center = roll_mean(roll_assets["x"], nwin, roll_assets["ind"], "center")

    left_nan_window = nwin * -1 + 1
    right_nan_window = nwin - 1

    left_grab_reals = [np.isnan(x) for x in left[: len(roll_assets["x"]) + left_nan_window]]
    right_grab_reals = [np.isnan(x) for x in right[right_nan_window:]]

    center_grab_reals_l = len(roll_assets["x"]) + int(floor(left_nan_window / 2.0))
    center_grab_reals_r = int(ceil(right_nan_window / 2.0))
    center_grab_reals = [np.isnan(x) for x in center[center_grab_reals_r:center_grab_reals_l]]

    assert all(not x for x in left_grab_reals)
    assert all(not x for x in right_grab_reals)
    assert all(not x for x in center_grab_reals)
