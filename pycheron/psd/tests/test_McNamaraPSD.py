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
from obspy.core.utcdatetime import UTCDateTime
import pytest

from pycheron.psd.McNamaraPSD import McNamaraPSD


@pytest.fixture(scope="module")
def trace(file_assets):
    return file_assets["test_traces"][0]


@pytest.mark.parametrize("binned", [True, False])
def test_McNamaraPSD_returns_correct_types(trace, binned):
    results = McNamaraPSD(trace)
    assert len(results) == 5
    # frequencies
    assert isinstance(results[0], np.ndarray)
    assert isinstance(results[0][0], np.float64)
    # spectrum
    assert isinstance(results[1], np.ndarray)
    assert isinstance(results[1][0], np.float64)
    # snlcq
    assert isinstance(results[2], list)
    assert isinstance(results[2][0], str)
    # Starttime string
    assert isinstance(results[3], str)
    # Endtime string
    assert isinstance(results[4], str)


def test_McNamaraPSD_returns_matching_start_end_times(trace):
    _, _, _, start_time, end_time = McNamaraPSD(trace)
    assert UTCDateTime(start_time) == trace.stats.starttime
    assert UTCDateTime(end_time) == trace.stats.endtime


def test_McNamaraPSD_returns_correct_snlcq_string_in_array_with_quality(trace):
    _, _, snlcq_string, _, _ = McNamaraPSD(trace)
    tr = trace
    expected = "{0}.{1}.{2}.{3}.{4}".format(
        str(tr.stats.network),
        str(tr.stats.station),
        str(tr.stats.location),
        str(tr.stats.channel),
        str(tr.stats.mseed.dataquality),
    )
    assert snlcq_string[0] == expected


def test_McNamaraPSD_returns_correct_snlcq_string_in_array_with_out_quality(trace):
    # File types other than MSEED may not have quality flag.
    # For simplicity, just remove the MSEED info from the trace instead of
    # loading another trace to test.
    trace.stats.pop("mseed")

    _, _, snlcq_string, _, _ = McNamaraPSD(trace)
    tr = trace
    expected = "{0}.{1}.{2}.{3}".format(
        str(tr.stats.network),
        str(tr.stats.station),
        str(tr.stats.location),
        str(tr.stats.channel),
    )
    assert snlcq_string[0] == expected
