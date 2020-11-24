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
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from pycheron.plotting.rankingPlot import rankPlot
from pycheron.db.sqllite_db import Database
from pycheron.test.create_test_db import create_test_db


@pytest.mark.parametrize("model", ["nlnm", "nnm"])
@pytest.mark.parametrize("f_name", [None, "ftest"])
@pytest.mark.parametrize("chan", ["BHN", "BHZ", "BHE"])
def test_st_rankPlot(f_name, model, chan):
    client = Client("IRIS")
    starttime = UTCDateTime("2012-12-12T00:00:00.000")
    endtime = UTCDateTime("2012-12-13T00:00:00.000")
    st = client.get_waveforms("AK", "GHO", "", chan, starttime, endtime)

    rval = rankPlot(st, model=model, f_name=f_name)
    assert rval is None

    if f_name is None:
        assert os.path.exists("f_name" + chan + "_1_AK.png")
        assert os.path.exists("f_name" + chan + "_6.5_AK.png")
        assert os.path.exists("f_name" + chan + "_30_AK.png")
        assert os.path.exists("f_name" + chan + "_100_AK.png")
        os.remove("f_name" + chan + "_1_AK.png")
        os.remove("f_name" + chan + "_6.5_AK.png")
        os.remove("f_name" + chan + "_30_AK.png")
        os.remove("f_name" + chan + "_100_AK.png")
    else:
        assert os.path.exists(f_name + chan + "_1_AK.png")
        assert os.path.exists(f_name + chan + "_6.5_AK.png")
        assert os.path.exists(f_name + chan + "_30_AK.png")
        assert os.path.exists(f_name + chan + "_100_AK.png")
        os.remove(f_name + chan + "_1_AK.png")
        os.remove(f_name + chan + "_6.5_AK.png")
        os.remove(f_name + chan + "_30_AK.png")
        os.remove(f_name + chan + "_100_AK.png")
    plt.close("all")


@pytest.mark.parametrize("model", ["nlnm", "nnm"])
@pytest.mark.parametrize("f_name", [None, "ftest"])
@pytest.mark.parametrize("chan", ["BHN", "BHZ", "BHE"])
def test_db_rankPlot(f_name, model, chan):
    create_test_db(sql_file="test_data/multi_chan_AK.sql")
    db_path = "test_data/pycheron_test_db.db"
    db = Database(db_path, "WFpych", True)

    rval = rankPlot(
        db, model=model, f_name=f_name, network="AK", station="GHO", channel=chan
    )
    assert rval is None

    if f_name is None:
        assert os.path.exists("f_name" + chan + "_1_AK.png")
        assert os.path.exists("f_name" + chan + "_6.5_AK.png")
        assert os.path.exists("f_name" + chan + "_30_AK.png")
        assert os.path.exists("f_name" + chan + "_100_AK.png")
        os.remove("f_name" + chan + "_1_AK.png")
        os.remove("f_name" + chan + "_6.5_AK.png")
        os.remove("f_name" + chan + "_30_AK.png")
        os.remove("f_name" + chan + "_100_AK.png")
    else:
        assert os.path.exists(f_name + chan + "_1_AK.png")
        assert os.path.exists(f_name + chan + "_6.5_AK.png")
        assert os.path.exists(f_name + chan + "_30_AK.png")
        assert os.path.exists(f_name + chan + "_100_AK.png")
        os.remove(f_name + chan + "_1_AK.png")
        os.remove(f_name + chan + "_6.5_AK.png")
        os.remove(f_name + chan + "_30_AK.png")
        os.remove(f_name + chan + "_100_AK.png")
    plt.close("all")


@pytest.mark.parametrize("model", ["nnm"])
@pytest.mark.parametrize("f_name", [None])
@pytest.mark.parametrize("chan", ["BHE"])
def test_st_rankPlot_with_db(f_name, model, chan):
    client = Client("IRIS")
    starttime = UTCDateTime("2012-12-12T00:00:00.000")
    endtime = UTCDateTime("2012-12-13T00:00:00.000")
    st = client.get_waveforms("AK", "GHO", "", chan, starttime, endtime)

    create_test_db(sql_file="test_data/multi_chan_AK.sql")
    db_path = "test_data/pycheron_test_db.db"
    db = Database(db_path, "WFpych", True)

    rval = rankPlot(
        st,
        model=model,
        f_name=f_name,
        network="AK",
        station="GHO",
        channel=chan,
        database=db,
    )
    assert rval is None
    assert os.path.exists("f_name" + chan + "_1_AK.png")
    assert os.path.exists("f_name" + chan + "_6.5_AK.png")
    assert os.path.exists("f_name" + chan + "_30_AK.png")
    assert os.path.exists("f_name" + chan + "_100_AK.png")
    os.remove("f_name" + chan + "_1_AK.png")
    os.remove("f_name" + chan + "_6.5_AK.png")
    os.remove("f_name" + chan + "_30_AK.png")
    os.remove("f_name" + chan + "_100_AK.png")
    plt.close("all")
