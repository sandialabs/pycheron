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
from obspy.core.stream import Stream
from obspy.core.trace import Trace


from pycheron.psd.noise.deadChannel import DDT, isDC


def _helper_construct_trace_with_gap(first_trace):
    old_stats = first_trace.stats
    second_trace_data = deepcopy(first_trace.data)
    second_trace_header = {}
    second_trace_header["sampling_rate"] = old_stats.sampling_rate
    second_trace_header["calib"] = old_stats.calib
    second_trace_header["npts"] = old_stats.npts
    second_trace_header["delta"] = old_stats.delta

    time_diff_original_trace = first_trace.stats.endtime - first_trace.stats.starttime
    second_trace_header["starttime"] = first_trace.stats.endtime + 10
    second_trace_header["endtime"] = (
        second_trace_header["starttime"] + time_diff_original_trace
    )
    second_trace = Trace(header=second_trace_header, data=second_trace_data)
    second_trace.id = first_trace.id
    first_trace += second_trace
    return first_trace


@pytest.mark.parametrize("dm", [True, False])
@pytest.mark.parametrize("dt", [True, False])
@pytest.mark.parametrize("tap", [0, 0.1, 0.5])
def test_DDT_returns_trace_obj_of_same_length_and_type(dm, dt, tap, file_assets):
    tr = file_assets["test_traces"][3].copy()
    result = DDT(tr, detrend=dt, demean=dm, taper=tap)
    assert len(result.data) == len(tr.data)
    original_type = type(tr.data[0])
    assert isinstance(result, Trace)
    assert isinstance(result.data[0], original_type)


@pytest.mark.parametrize("dm", [True, False])
@pytest.mark.parametrize("dt", [True, False])
@pytest.mark.parametrize("tap", [0, 0.1, 0.5])
def test_DDT_returns_trace_obj_for_single_trace_streams(dm, dt, tap, file_assets):
    st = file_assets["test_streams"][3].copy()
    result = DDT(st, detrend=dt, demean=dm, taper=tap)
    assert len(result.data) == len(st[0].data)
    assert isinstance(result, Trace)


@pytest.mark.parametrize("dm", [True, False])
@pytest.mark.parametrize("dt", [True, False])
@pytest.mark.parametrize("tap", [0, 0.1, 0.5])
def test_DDT_returns_trace_obj_for_traces_with_gaps(file_assets, dm, dt, tap):
    # Trace must be added to another trace with a gap in time between the two to
    # hit this case.
    old_trace = file_assets["test_traces"][0].copy()
    trace = _helper_construct_trace_with_gap(old_trace)
    result = DDT(trace, detrend=dt, demean=dm, taper=tap)
    assert isinstance(result, Trace)


@pytest.mark.parametrize("dm", [True, False])
@pytest.mark.parametrize("dt", [True, False])
@pytest.mark.parametrize("tap", [0, 0.1, 0.5])
def test_DDT_returns_trace_obj_for_single_trace_stream_with_gaps(
    file_assets, dm, dt, tap
):
    # Trace must be added to another trace with a gap in time between the two to
    # hit this case.
    old_trace = file_assets["test_traces"][3].copy()
    trace = _helper_construct_trace_with_gap(old_trace)
    stream = Stream([trace])
    result = DDT(stream, detrend=dt, demean=dm, taper=tap)
    assert isinstance(result, Trace)


@pytest.mark.parametrize("dm", [True, False])
@pytest.mark.parametrize("dt", [True, False])
@pytest.mark.parametrize("tap", [0, 0.1, 0.5])
def test_isDC_returns_false_for_live_channel_trace(dm, dt, tap, file_assets):
    tr = file_assets["test_traces"][0].copy()
    result = isDC(tr, detrend=dt, demean=dm, taper=tap)
    assert result == False


@pytest.mark.parametrize("dm", [True, False])
@pytest.mark.parametrize("dt", [True, False])
@pytest.mark.parametrize("tap", [0, 0.1, 0.5])
def test_isDC_returns_true_for_trace_of_all_zeros_trace(dm, dt, tap, file_assets):
    tr = file_assets["test_traces"][3].copy()
    tr.data = np.zeros(len(tr.data))
    result = isDC(tr, detrend=dt, demean=dm, taper=tap)
    assert result == True


@pytest.mark.parametrize("dm", [True, False])
@pytest.mark.parametrize("dt", [True, False])
@pytest.mark.parametrize("tap", [0, 0.1, 0.5])
def test_isDC_returns_true_for_trace_of_all_zeros_stream(dm, dt, tap, file_assets):
    st = file_assets["test_streams"][3].copy()
    st[0].data = np.zeros(len(st[0].data))
    result = isDC(st, detrend=dt, demean=dm, taper=tap)
    assert result == True


@pytest.mark.parametrize("dm", [True, False])
@pytest.mark.parametrize("dt", [True, False])
@pytest.mark.parametrize("tap", [0, 0.1, 0.5])
def test_isDC_returns_false_for_live_channel_stream(dm, dt, tap, file_assets):
    st = file_assets["test_streams"][3].copy()
    result = isDC(st, detrend=dt, demean=dm, taper=tap)
    assert result == False
