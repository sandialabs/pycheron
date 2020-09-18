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
from pycheron.rollseis.roll_median import _medianFilter, roll_median


@pytest.mark.parametrize("x", [np.array([1, 2, 3]), []])
def test_nwin_gt_x_roll_median(roll_assets, x):
    median_vals = roll_median(x, roll_assets["nwin"], roll_assets["ind"])
    assert median_vals is None


@pytest.mark.parametrize("ind", list(range(-10, 1)))
def test_ind_lt_1_roll_median(roll_assets, ind):
    median_vals = roll_median(roll_assets["x"], roll_assets["nwin"], ind)
    assert median_vals is None


@pytest.mark.parametrize("nwin", [6, 7])
@pytest.mark.parametrize("ind", [1, 2, 3, 4, 5])
def test_nans_median_roll(roll_assets, nwin, ind):
    median_vals = roll_median(roll_assets["x"], nwin, ind)
    # for a median_vals window: the first nwin-1 / 2 and last nwin-1 / 2 values should be NaN
    # (nwin increased by 1 in order to guarantee odd window size if necessary)
    left_nan_window = int(floor((nwin * -1 + 1) / 2.0))
    right_nan_window = int(ceil((nwin - 1) / 2.0))

    median_vals_grab_nans_l = [np.isnan(x) for x in median_vals[left_nan_window:]]
    median_vals_grab_nans_r = [np.isnan(x) for x in median_vals[:right_nan_window]]

    assert len(median_vals_grab_nans_l) == len(median_vals_grab_nans_r)

    assert all(x for x in median_vals_grab_nans_l)
    assert all(x for x in median_vals_grab_nans_r)


@pytest.mark.parametrize("nwin", [6, 7])
def test_reals_median_roll(roll_assets, nwin):
    median_vals = roll_median(roll_assets["x"], nwin, roll_assets["ind"])
    left_nan_window = int(floor((nwin * -1 + 1) / 2.0))
    right_nan_window = int(ceil((nwin - 1) / 2.0))

    median_vals_grab_reals_l = len(roll_assets["x"]) + int(floor(left_nan_window))
    median_vals_grab_reals_r = int(ceil(right_nan_window))
    median_vals_grab_reals = [
        np.isnan(x)
        for x in median_vals[median_vals_grab_reals_r:median_vals_grab_reals_l]
    ]

    assert all(not x for x in median_vals_grab_reals)
