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

from pycheron.sigpro.unHistogram import unHistogram


def is_sorted_increasing(list_to_check):
    return all(
        list_to_check[i] <= list_to_check[i + 1] for i in range(len(list_to_check) - 1)
    )


def test_length_four_increasing_array_in_unhistogram():
    results = unHistogram([1, 2, 3, 4])
    assert all(isinstance(x, float64) for x in results)
    assert is_sorted_increasing(results)
    assert len(results) == 10


def test_length_five_increasing_array_in_unhistogram():
    results = unHistogram([1, 2, 3, 4, 5])
    assert all(isinstance(x, float64) for x in results)
    assert is_sorted_increasing(results)
    assert len(results) == 470


def test_length_seven_increasing_array_in_unhistogram():
    results = unHistogram([1, 2, 3, 4, 5, 6, 7])
    assert all(isinstance(x, float64) for x in results)
    assert is_sorted_increasing(results)
    assert len(results) == 470
