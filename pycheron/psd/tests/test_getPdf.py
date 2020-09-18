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

from pycheron.psd.getPDF import getPDF
from pycheron.psd.noise.getNoise import getNoise
from pycheron.psd.psdList import psdList


@pytest.fixture(scope="module")
def noise_matrix(file_assets):
    test_stream = file_assets["test_streams"][2]
    # Setup code from example in getPDF.py
    # creating psds
    psd = psdList(test_stream)
    # getting noise, frequency
    f, n, psd = getNoise(psd)
    return n[0]


def test_getPDF_returns_numpy_array_of_arrays_of_floats(noise_matrix):
    result = getPDF(noise_matrix)
    assert isinstance(result, np.ndarray)
    assert isinstance(result[0], np.ndarray)
    assert isinstance(result[0][0], np.float64)


def test_getPDF_returns_array_of_length_diff_lo_hi_over_binSize(noise_matrix):
    result = getPDF(noise_matrix, lo=20, hi=40)
    assert len(result) == 20

    result = getPDF(noise_matrix, lo=-100, hi=-20)
    assert len(result) == 80

    result = getPDF(noise_matrix, lo=0, hi=0)
    assert len(result) == 0

    result = getPDF(noise_matrix, lo=20, hi=40, binSize=2)
    assert len(result) == 10
