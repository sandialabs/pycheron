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

from pycheron.util.getUpDownTimes import getUpDownTimes


def test_getUpDownTimes(file_assets):
    # TODO: The block to check delta value is never caught (delta is absolute value)
    # can this be removed?

    # uid > 1, with 1, 2, >2 traces
    up_down_times_multi = getUpDownTimes(file_assets["multi_stream"])
    assert len(up_down_times_multi) > 0
    assert all(isinstance(x, dict) for x in up_down_times_multi)
    assert all(bool(x) for x in up_down_times_multi)

    up_down_times_two = getUpDownTimes(file_assets["7A_CABN_stream"] + file_assets["BGU_HHN_stream"])
    assert len(up_down_times_two) > 0
    assert all(isinstance(x, dict) for x in up_down_times_two)
    assert all(bool(x) for x in up_down_times_two)

    up_down_times_7A = getUpDownTimes(file_assets["7A_CABN_stream"])
    assert len(up_down_times_7A) > 0
    assert all(isinstance(x, dict) for x in up_down_times_7A)
    assert all(bool(x) for x in up_down_times_7A)


def test_uid_eq1_getUpDownTimes(file_assets):
    # uid == 1, with 1, 2, >2 traces
    up_down_times_multi = getUpDownTimes(
        file_assets["7A_CABN_stream"] + file_assets["7A_CABN_stream"] + file_assets["7A_CABN_stream"]
    )
    assert len(up_down_times_multi) > 0
    assert all(isinstance(x, dict) for x in up_down_times_multi)
    assert all(bool(x) for x in up_down_times_multi)

    up_down_times_two = getUpDownTimes(file_assets["7A_CABN_stream"] + file_assets["7A_CABN_stream"])
    assert len(up_down_times_two) > 0
    assert all(isinstance(x, dict) for x in up_down_times_two)
    assert all(bool(x) for x in up_down_times_two)

    up_down_times_7A = getUpDownTimes(file_assets["7A_CABN_stream"])
    assert len(up_down_times_7A) > 0
    assert all(isinstance(x, dict) for x in up_down_times_7A)
    assert all(bool(x) for x in up_down_times_7A)


def test_1_stream_min_signal_gt_signal_duration_getUpDownTimes(file_assets):
    up_down_times_7A = getUpDownTimes(file_assets["7A_CABN_stream"], min_signal=999999999)
    assert up_down_times_7A == []


def test_3_streams_min_signal_gt_signal_duration_getUpDownTimes(file_assets):
    up_down_times_multi = getUpDownTimes(file_assets["multi_stream"], min_signal=999999999)
    assert up_down_times_multi == []


def test_delta_lt_min_gap_getUpDownTimes(file_assets):
    # uid > 1, with 2, >2 traces
    up_down_times_multi = getUpDownTimes(file_assets["multi_stream"], min_gap=999999999)
    for vals in up_down_times_multi:
        assert vals["start_time"] == "NaN"
        assert vals["end_time"] == "NaN"
        assert vals["channel_up_time"] == 100

    up_down_times_two = getUpDownTimes(file_assets["7A_CABN_stream"] + file_assets["BGU_HHN_stream"], min_gap=999999999)
    for vals in up_down_times_two:
        assert vals["start_time"] == "NaN"
        assert vals["end_time"] == "NaN"
        assert vals["channel_up_time"] == 100


def test_uid_eq_1_delta_lt_min_gap_getUpDownTimes(file_assets):
    # uid == 1, with 2, >2 traces
    up_down_times_multi = getUpDownTimes(
        file_assets["7A_CABN_stream"] + file_assets["7A_CABN_stream"] + file_assets["7A_CABN_stream"],
        min_gap=999999999,
    )
    for vals in up_down_times_multi:
        assert vals["channel_up_time"] == 100

    for vals in up_down_times_multi[1:]:
        assert vals["start_time"] == "NaN"
        assert vals["end_time"] == "NaN"

    up_down_times_two = getUpDownTimes(file_assets["7A_CABN_stream"] + file_assets["7A_CABN_stream"], min_gap=999999999)
    for vals in up_down_times_two:
        assert vals["channel_up_time"] == 100

    for vals in up_down_times_two[1:]:
        assert vals["start_time"] == "NaN"
        assert vals["end_time"] == "NaN"
