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
from pycheron.util.logger import Logger


@pytest.fixture(scope="module")
def cross_spectrum_func_w_data(file_assets):
    test_trace_1 = file_assets["test_traces"][0]
    test_trace_2 = file_assets["test_traces"][1]
    sampling_rate = int(round(test_trace_1.stats.sampling_rate))

    def cross_spectrum_call(
        spans=np.array([3, 5, 7, 9]),
        taper=0.1,
        pad=0,
        demean=False,
        detrend=True,
        logger=None,
    ):
        return crossSpectrum(
            test_trace_1,
            test_trace_2,
            sampling_rate,
            spans=spans,
            taper=taper,
            pad=pad,
            demean=demean,
            detrend=detrend,
            logger=logger,
        )

    return cross_spectrum_call


def assert_is_non_empty_numpy_array_with_type(test_value, type):
    # Helper Function for asserting numpy array types
    assert isinstance(test_value, np.ndarray)
    assert isinstance(test_value[0], type)


def assert_return_of_cross_spectrum_has_correct_types(return_arrays):
    # Helper function for casserting cross spectrum function returns
    # correct length array with correct types.
    assert len(return_arrays) == 8

    spectral_frequencies = return_arrays[0]
    assert_is_non_empty_numpy_array_with_type(spectral_frequencies, np.float64)

    trace_1_spectral_amplitudes = return_arrays[1]
    assert_is_non_empty_numpy_array_with_type(trace_1_spectral_amplitudes, np.float64)

    trace_2_spectral_amplitudes = return_arrays[2]
    assert_is_non_empty_numpy_array_with_type(trace_2_spectral_amplitudes, np.float64)

    magnitude_squared_coherence = return_arrays[3]
    assert_is_non_empty_numpy_array_with_type(magnitude_squared_coherence, np.float64)

    cross_spectral_phase = return_arrays[4]
    assert_is_non_empty_numpy_array_with_type(cross_spectral_phase, np.float64)

    trace_1_periodogram = return_arrays[5]
    assert_is_non_empty_numpy_array_with_type(trace_1_periodogram, np.complex128)

    trace_2_periodogram = return_arrays[6]
    assert_is_non_empty_numpy_array_with_type(trace_2_periodogram, np.complex128)

    cross_periodogram = return_arrays[7]
    assert_is_non_empty_numpy_array_with_type(cross_periodogram, np.complex128)


def test_crossSpectrum_default_params_returns_correct_types(cross_spectrum_func_w_data):
    results = cross_spectrum_func_w_data()
    assert_return_of_cross_spectrum_has_correct_types(results)


@pytest.mark.parametrize("demean", [True, False])
@pytest.mark.parametrize("detrend", [True, False])
def test_demean_detrend_return_correct_types(demean, detrend, cross_spectrum_func_w_data):
    results = cross_spectrum_func_w_data(demean=demean, detrend=detrend)
    assert_return_of_cross_spectrum_has_correct_types(results)


def test_zero_taper_returns_correct_types(cross_spectrum_func_w_data):
    results = cross_spectrum_func_w_data(taper=0)
    assert_return_of_cross_spectrum_has_correct_types(results)


def test_pad_greater_than_zero_returns_correct_types(cross_spectrum_func_w_data):
    results = cross_spectrum_func_w_data(pad=1)
    assert_return_of_cross_spectrum_has_correct_types(results)


def test_crossSpectrum_with_no_span_returns_None(cross_spectrum_func_w_data):
    # This should raise an exception in the future. Currently logs error and
    # returns None.
    result = cross_spectrum_func_w_data(spans=None, logger=Logger(None))
    assert result is None
