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

from pycheron.sigpro.cosine_bell import cosine_bell


def test_cosine_bell_returns_ndarray_of_floats(file_assets):
    result = cosine_bell(file_assets["test_trace"].data)
    assert isinstance(result, np.ndarray)
    assert isinstance(result[0], np.float64)


def test_cosine_bell_with_different_taper_returns_different_values(file_assets):
    result_one_tenth = cosine_bell(file_assets["test_trace"].data, 0.1)
    result_two_tenth = cosine_bell(file_assets["test_trace"].data, 0.2)
    assert (result_one_tenth != result_two_tenth).any()


def test_taper_less_than_1_or_greater_than_one_half_returns_none(file_assets):
    # Currently, when an invalid value is passed into the taper parameter,
    # the function returns None after logging an error. This should probably be
    # changed to raise a ValueError.
    result_negative = cosine_bell(file_assets["test_trace"].data, -1)
    result_greater_than_half = cosine_bell(file_assets["test_trace"].data, 0.51)
    assert result_negative is None
    assert result_greater_than_half is None
