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

from copy import deepcopy
import numpy as np
import pytest

from pycheron.psd.psdList import psdList
from pycheron.psd.noise.getNoise import getNoise, _pop_indices
from pycheron.util.logger import Logger


@pytest.fixture(scope="module")
def psd_list_data(file_assets):
    st = file_assets["test_streams"][3]
    return psdList(st)


def test_get_noise_noise_matricies_are_correct_types_valid_trace(psd_list_data):
    f, n, psd = getNoise(psd_list_data)
    assert isinstance(f, list)
    assert isinstance(f[0], np.ndarray)
    assert isinstance(f[0][0], float)
    assert isinstance(n, list)
    assert isinstance(n[0], np.ndarray)
    assert isinstance(n[0][0], np.ndarray)
    assert isinstance(n[0][0][0], float)


def test_get_noise_noise_matricies_return_correct_values(psd_list_data):
    f, n, psd = getNoise(psd_list_data)
    assert f[0][0] == pytest.approx(0.005255602595335718, 0.0001)
    assert f[0][-1] == pytest.approx(9.870149282610823, 0.0001)
    assert n[0][0][0] == pytest.approx(-126.33487879809951, 0.001)
    assert n[0][0][-1] == pytest.approx(-82.64411867131393, 0.001)
    assert psd[0][0][0][0] == pytest.approx(0.005255602595335718, 0.0001)
    assert psd[0][0][0][-1] == pytest.approx(9.870149282610823, 0.0001)


def test_get_noise_noise_matricies_are_correct_types_trace_w_out_quality(psd_list_data):
    test_data = deepcopy(psd_list_data)
    # get rid of quality metadata from MSEED
    snlcq_array = str(test_data[0][0][2][0]).split(".")
    snlcq_array.pop(4)
    snlcq_string = ".".join(snlcq_array)
    test_data[0][0][2][0] = snlcq_string
    f, n, psd = getNoise(test_data)
    assert isinstance(f, list)
    assert isinstance(f[0], np.ndarray)
    assert isinstance(f[0][0], float)
    assert isinstance(n, list)
    assert isinstance(n[0], np.ndarray)
    assert isinstance(n[0][0], np.ndarray)
    assert isinstance(n[0][0][0], float)


def test_get_noise_pops_psdLists_without_valid_information(psd_list_data):
    test_data = deepcopy(psd_list_data)
    # get rid of quality metadata from MSEED
    snlcq_array = str(test_data[0][0][2][0]).split(".")
    snlcq_array[1] = "nonsense station"
    snlcq_string = ".".join(snlcq_array)
    test_data[0][0][2][0] = snlcq_string
    f, n, psd = getNoise(test_data, logger=Logger(None))
    assert isinstance(f, list)
    assert len(f) == 0
    assert isinstance(n, list)
    assert len(n) == 0
    assert isinstance(psd, list)
    assert len(psd) == 0


def test_get_noise_removes_non_list_types_from_psdlist(psd_list_data):
    # TODO: Figure out if we really need this functionality. If we do need it,
    # improve it to handle strings, dictionaries, and other classes with a
    # length attribute. This function tries to call the len function
    # on every element of the psdList. If the len call throws a TypeError, then
    # it pops that element. Since strings have a len attribute, they will pass
    # this check, but error later on. Same for dictionaries, etc... Should
    # probably be checking to see if elements are a list or a numpy array.
    test_data = deepcopy(psd_list_data)
    test_data.append(None)
    test_data.append(1)
    test_data.append(True)
    f, n, psd = getNoise(test_data)
    assert isinstance(f, list)
    assert len(f) == 1
    assert isinstance(n, list)
    assert isinstance(psd, list)
    assert len(psd) == 1


def test_pop_indicies_pos_supplied_indicies():
    # This function doesn't work with negative indicies.
    indicies_to_pop = [0, 2, 5]
    values = [i for i in range(6)]
    result = _pop_indices(indicies_to_pop, values)
    assert result == [1, 3, 4]
