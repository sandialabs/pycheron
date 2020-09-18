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

from numpy import float64
import pytest

from pycheron.sigpro.demeanTaper import cosine_taper


def test_quarter_cosine_taper_returns_np_float_array_of_length_npts():
    results = cosine_taper(npts=40, halfcosine=False)
    assert all(isinstance(x, float64) for x in results)
    assert len(results) == 40


def test_quarter_cosine_SAC_taper_returns_np_float_array_of_length_npts():
    results = cosine_taper(npts=40, halfcosine=False, sactaper=True)
    assert all(isinstance(x, float64) for x in results)
    assert len(results) == 40


def test_half_cosine_taper_returns_np_float_array_of_length_npts():
    results = cosine_taper(npts=40, halfcosine=True)
    assert all(isinstance(x, float64) for x in results)
    assert len(results) == 40


def test_half_cosine_SAC_taper_returns_np_float_array_of_length_npts():
    results = cosine_taper(npts=40, halfcosine=True)
    assert all(isinstance(x, float64) for x in results)
    assert len(results) == 40


def test_percentage_equal_to_1_or_0_taper_maintains_above_functionalities():
    results_p_1 = cosine_taper(npts=40, p=1.0)
    assert all(isinstance(x, float64) for x in results_p_1)
    assert len(results_p_1) == 40

    results_p_2 = cosine_taper(npts=40, p=0.0)
    assert all(isinstance(x, float64) for x in results_p_2)
    assert len(results_p_2) == 40


def test_idx1_equals_idx2_maintains_above_functionalities():
    results = cosine_taper(npts=10)
    assert all(isinstance(x, float64) for x in results)
    assert len(results) == 10


def test_npts_equal_to_zero_returns_length_zero_array():
    results = cosine_taper(npts=0)
    assert len(results) == 0


def test_p_less_than_0_or_greater_than_1_raises_value_error():
    with pytest.raises(ValueError):
        cosine_taper(npts=1, p=1.01)
        cosine_taper(npts=1, p=-0.01)
