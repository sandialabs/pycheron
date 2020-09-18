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

from pycheron.psd.noise.noiseModel import noiseModel
from pycheron.psd.noise.getNoise import getNoise


@pytest.fixture(scope="module")
def frequencies(BHZ_channel_stream_and_psd_list):
    psds = BHZ_channel_stream_and_psd_list[1]
    f, n, psd = getNoise(psds)
    freq = f[0]
    return freq


@pytest.fixture(scope="module")
def low_and_high_noise_models(frequencies):
    return noiseModel(frequencies)


def noiseModel_returns_two_numpy_arrays(models):
    assert len(models) == 2
    assert all([isinstance(model, np.ndarray) for model in models])


def test_noiseModel_returns_two_numpy_arrays(low_and_high_noise_models):
    noiseModel_returns_two_numpy_arrays(low_and_high_noise_models)


def test_low_noise_is_numpy_floats(low_and_high_noise_models):
    low_noise = low_and_high_noise_models[0]
    assert all([isinstance(noise, np.float64) for noise in low_noise])


def test_low_noise_values(low_and_high_noise_models):
    low_noise = low_and_high_noise_models[0]
    assert low_noise[0] == pytest.approx(-185.68606263476232, 0.01)
    assert low_noise[-1] == pytest.approx(-167.96798578798467, 0.01)


def test_high_noise_is_numpy_floats(low_and_high_noise_models):
    low_noise = low_and_high_noise_models[1]
    assert all([isinstance(noise, np.float64) for noise in low_noise])


def test_high_noise_values(low_and_high_noise_models):
    high_noise = low_and_high_noise_models[1]
    assert high_noise[0] == pytest.approx(-128.70343140946508, 0.01)
    assert high_noise[-1] == pytest.approx(-91.59780228245113, 0.01)


def test_noiseModel_returns_numpy_array_if_passed_list(frequencies):
    frequencies = frequencies.tolist()
    low_and_high_models = noiseModel(frequencies)
    noiseModel_returns_two_numpy_arrays(low_and_high_models)
