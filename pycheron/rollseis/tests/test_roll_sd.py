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
from math import ceil, floor
from pycheron.rollseis.roll_sd import _sdFilter, roll_sd


@pytest.mark.parametrize("nwin", [6, 7])
@pytest.mark.parametrize("ind", [1, 2, 3, 4, 5])
def test_nans_sd_roll(roll_assets, nwin, ind):
    # TODO: for refactor, ensure "left", "right", and "center" are the only values accepted
    # currently, "right" will produce the same result as a random string
    #
    # For refactor:
    # There needs to be a check for when ind<0 in the original code
    # This will result in an infinite loop
    # There needs to be a check for when nwin > len(x), this breaks the function
    left = roll_sd(roll_assets["x"], nwin, ind, "left")
    right = roll_sd(roll_assets["x"], nwin, ind, "right")
    center = roll_sd(roll_assets["x"], nwin, ind, "center")

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

    center_grab_nans_l = [
        np.isnan(x) for x in center[int(floor(left_nan_window / 2.0)) :]
    ]
    center_grab_nans_r = [
        np.isnan(x) for x in center[: int(ceil(right_nan_window / 2.0))]
    ]

    assert len(center_grab_nans_l) == len(center_grab_nans_r)

    assert all(x for x in left_grab_nans)
    assert all(x for x in right_grab_nans)
    assert all(x for x in center_grab_nans_l)
    assert all(x for x in center_grab_nans_r)


@pytest.mark.parametrize("nwin", [6, 7])
def test_reals_sd_roll(roll_assets, nwin):
    # if ind == 1, the only NaN values expected should be:
    # for a left window: the last nwin-1 values
    # for a right window: the first nwin-1 values
    # for a center window: the first nwin-1 / 2 and last nwin-1 / 2 values should be NaN
    # (nwin increased by 1 in order to guarantee odd window size if necessary)
    # all others should real numbers
    left = roll_sd(roll_assets["x"], nwin, roll_assets["ind"], "left")
    right = roll_sd(roll_assets["x"], nwin, roll_assets["ind"], "right")
    center = roll_sd(roll_assets["x"], nwin, roll_assets["ind"], "center")

    left_nan_window = nwin * -1 + 1
    right_nan_window = nwin - 1

    left_grab_reals = [
        np.isnan(x) for x in left[: len(roll_assets["x"]) + left_nan_window]
    ]
    right_grab_reals = [np.isnan(x) for x in right[right_nan_window:]]

    center_grab_reals_l = len(roll_assets["x"]) + int(floor(left_nan_window / 2.0))
    center_grab_reals_r = int(ceil(right_nan_window / 2.0))
    center_grab_reals = [
        np.isnan(x) for x in center[center_grab_reals_r:center_grab_reals_l]
    ]

    assert all(not x for x in left_grab_reals)
    assert all(not x for x in right_grab_reals)
    assert all(not x for x in center_grab_reals)
