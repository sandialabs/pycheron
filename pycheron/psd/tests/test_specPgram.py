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

from pycheron.psd.specPgram import specPgram


def assert_specPgram_types_are_correct(results):
    # Helper function to assert correct return types of specPgram.
    assert len(results) == 2
    for result in results:
        assert isinstance(result, np.ndarray)
        assert isinstance(result[0], np.float64)


@pytest.mark.parametrize("demean", [True, False])
@pytest.mark.parametrize("detrend", [True, False])
@pytest.mark.parametrize("fast", [True, False])
@pytest.mark.parametrize("taper", [0, 0.1, 0.5])
def test_specPgram_returns_correct_types(file_assets, demean, detrend, fast, taper):
    trace = file_assets["test_traces"][3]
    sampling_frequency = trace.stats.sampling_rate
    results = specPgram(
        trace,
        sampling_frequency,
        taper=taper,
        fast=fast,
        demean=demean,
        detrend=detrend,
    )
    assert_specPgram_types_are_correct(results)
