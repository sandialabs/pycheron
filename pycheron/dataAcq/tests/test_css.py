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

import os
import pytest
import numpy as np
from obspy.core.stream import Stream
from obspy.core.trace import Trace
from obspy.core.utcdatetime import UTCDateTime
from pisces.io.readwaveform import read_waveform
from pycheron.dataAcq.css import (
    read_wfdisc,
    wf_dir_path_create,
    get_wfdisc_stations,
    create_wfdisc_header,
    css2stream,
    one_stream_per_day,
)


@pytest.fixture
def wfdisc_assets():
    return {"dir": "test_data/", "individual": "dpg_jun_17_2010_tests.wfdisc"}


def test_read_wfdisc(wfdisc_assets):
    good_wfdisc = wfdisc_assets["dir"] + wfdisc_assets["individual"]
    columns = [
        "sta",
        "chan",
        "time",
        "wfid",
        "chanid",
        "jdate",
        "endtime",
        "nsamp",
        "samprate",
        "calib",
        "calper",
        "instype",
        "segtype",
        "datatype",
        "clip",
        "dir",
        "dfile",
        "foff",
        "commid",
        "lddate",
    ]

    valid_df = read_wfdisc(good_wfdisc)
    assert valid_df.columns.tolist() == columns
    assert valid_df.empty is False
    assert valid_df.dropna().empty is False


def test_wf_dir_path_create(wfdisc_assets):
    valid_path = wf_dir_path_create(wfdisc_assets["dir"])
    assert valid_path.split(".")[-1] == "wfdisc"
    assert valid_path.split("/")[0] == wfdisc_assets["dir"].split("/")[0]


def test_get_wfdisc_stations(wfdisc_assets):
    wfdisc_dir = wfdisc_assets["dir"]
    unique_stations = get_wfdisc_stations(wfdisc_dir)

    assert isinstance(unique_stations, np.ndarray)
    assert len(unique_stations) > 0


def test_create_wfdisc_header(wfdisc_assets):
    good_wfdisc = wfdisc_assets["dir"] + wfdisc_assets["individual"]
    expected_dict_keys = [
        "starttime",
        "npts",
        "calib",
        "sampling_rate",
        "endtime",
        "channel",
        "network",
        "station",
    ]
    valid_df = read_wfdisc(good_wfdisc)

    valid_header = create_wfdisc_header(valid_df.loc[0], "TEST_NETWORK")
    assert isinstance(valid_header, dict)

    sorted_header_keys = list(valid_header.keys()).sort()
    sorted_expected_dict_keys = expected_dict_keys.sort()
    assert sorted_header_keys == sorted_expected_dict_keys
    assert isinstance(valid_header["starttime"], UTCDateTime)
    assert isinstance(valid_header["npts"], np.int64)
    assert isinstance(valid_header["calib"], np.float64)
    assert isinstance(valid_header["sampling_rate"], np.float64)
    assert isinstance(valid_header["endtime"], UTCDateTime)
    assert isinstance(valid_header["channel"], str)
    assert isinstance(valid_header["network"], str)
    assert isinstance(valid_header["station"], str)


@pytest.mark.parametrize(
    "station",
    np.array(
        [
            "DP91I1",
            "DP91I2",
            "DP91I3",
            "DP91I4",
            "DP91S1",
            "DP92I1",
            "DP92I2",
            "DP92I3",
            "DP92I4",
            "DP92S1",
        ],
        dtype=object,
    ),
)
def test_one_stream_per_day(wfdisc_assets, station):
    wf_file = wf_dir_path_create(wfdisc_assets["dir"])
    wfdisc = read_wfdisc(wf_file)
    # making traces
    traces = []
    starttimes_list = []
    # loop through column names
    for i, row in wfdisc.iterrows():
        if station == row["sta"]:
            header = create_wfdisc_header(row, "TEST_NETWORK")
            starttimes_list.append(UTCDateTime(row["time"]).isoformat()[0:10])
            wfile = os.path.join(wfdisc_assets["dir"], row["dfile"])

            # reading waveform
            data = read_waveform(wfile, row["datatype"], row["foff"], row["nsamp"])

            tr = Trace(data, header)
            traces.append(tr)
    one_per_day = one_stream_per_day(starttimes_list, traces)

    assert len(one_per_day) > 0
    assert all(isinstance(x, Stream) for x in one_per_day)


@pytest.mark.parametrize("byDay", [True, False])
@pytest.mark.parametrize(
    "station",
    np.array(
        [
            "DP91I1",
            "DP91I2",
            "DP91I3",
            "DP91I4",
            "DP91S1",
            "DP92I1",
            "DP92I2",
            "DP92I3",
            "DP92I4",
            "DP92S1",
        ],
        dtype=object,
    ),
)
def test_css2stream(wfdisc_assets, byDay, station):
    wfdisc_dir = wfdisc_assets["dir"]
    network = "TEST_NETWORK"

    valid_stream = css2stream(wfdisc_dir, network, station, byDay)

    if byDay is True:
        assert isinstance(valid_stream[0], Stream)
    else:
        assert isinstance(valid_stream, Stream)
