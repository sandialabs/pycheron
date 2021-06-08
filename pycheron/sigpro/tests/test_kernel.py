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

from pycheron.sigpro import kernel
from pycheron.util.logger import Logger


def test_one_dimensional_kernel():
    kernel_coefs = kernel.modified_daniell_kernel(1)
    assert kernel_coefs == [0.25, 0.5, 0.25]


def test_two_dimensional_kernel():
    kernel_coefs = kernel.modified_daniell_kernel(2)
    assert kernel_coefs == [0.125, 0.25, 0.25, 0.25, 0.125]


def test_multi_dimensional_span():
    # pytest 5.3.0 and higher support complex number tolerances in the
    # pytest.approx module. This should probably be set up as a +- small
    # tolerance test with actual values (see example in kernel.py for values)
    # instead of just testing if the number is complex.
    # see https://github.com/pytest-dev/pytest/issues/6057 for complex numbers
    spans = np.array([1, 2, 3, 4])
    kernel_coefs = kernel.modified_daniell_kernel(spans)

    assert len(kernel_coefs) == 21
    assert all(isinstance(x, complex) for x in kernel_coefs)


def test_kernapply_returns_x_when_dimension_0():
    result = kernel._kernapply(x=[1], k=None, m=0)
    assert result == [1]


def test_kernapply_returns_none_when_input_vector_less_than_2_times_dimension():
    # Currently, when an invalid length input vector is passed into the \
    # x parameter, the function returns None after logging an error.
    # This should probably be changed to raise a ValueError.
    result = kernel._kernapply(x=[1], k=None, m=5, logger=Logger(None))
    assert result is None
