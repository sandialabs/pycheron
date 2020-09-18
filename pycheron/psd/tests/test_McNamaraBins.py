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

from pycheron.psd.crossSpectrum import crossSpectrum
from pycheron.psd.McNamaraBins import McNamaraBins


@pytest.fixture(scope="module")
def frequency_and_spectral_values(file_assets):
    test_trace_1 = file_assets["test_traces"][0]
    test_trace_2 = file_assets["test_traces"][1]
    sampling_rate = int(round(test_trace_1.stats.sampling_rate))

    result = crossSpectrum(
        test_trace_1, test_trace_2, sampling_rate, spans=np.array([3, 5, 7, 9])
    )
    frequencies = result[0]
    spectral_amplitudes_trace_1 = result[1]

    return frequencies, spectral_amplitudes_trace_1


@pytest.mark.parametrize("align_freq", [0.004, 0.1, 1, 5])
def test_McNamaraBins_returns_two_numpy_arrays_of_floats(
    align_freq, frequency_and_spectral_values
):
    # Passing zero, negative numbers, or numbers larger than trace length in to
    # align_freq will cause uncaught exceptions. Probably should be changed in future.
    f = frequency_and_spectral_values[0]
    spec = frequency_and_spectral_values[1]

    results = McNamaraBins(f, spec, align_freq)
    assert isinstance(results[0], np.ndarray)
    assert isinstance(results[1], np.ndarray)
    assert isinstance(results[0][0], np.float64)
    assert isinstance(results[1][0], np.float64)


def test_McNamaraBins_align_freq_greater_than_hiFreq_returns_empty_array(
    frequency_and_spectral_values,
):
    f = frequency_and_spectral_values[0]
    spec = frequency_and_spectral_values[1]

    results = McNamaraBins(f, spec, 11)
    assert isinstance(results[0], np.ndarray)
    assert isinstance(results[1], np.ndarray)
    assert len(results[0]) == len(results[1]) == 0
