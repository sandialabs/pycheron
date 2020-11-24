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

import obspy
import os
from obspy.imaging.cm import pqlx
import pytest
import matplotlib.pyplot as plt
from pycheron.plotting.psdPlot import psdPlot
from pycheron.db.sqllite_db import Database
from pycheron.test.create_test_db import create_test_db


@pytest.mark.parametrize("style", ["psd", "pdf"])
@pytest.mark.parametrize("envelope_type", ["10_90", "05_95"])
@pytest.mark.parametrize("f_name", [None, "ftest"])
@pytest.mark.parametrize("single_type", ["5"])
@pytest.mark.parametrize("timespan", [0])
def test_psdPlot_7A(style, envelope_type, f_name, single_type, timespan):
    data = "test_data/7a_cabn_bhe.884965.tar.mseed"
    st = obspy.read(data)
    chan = "BHE"
    rval = psdPlot(
        st,
        style=style,
        f_name=f_name,
        showNoiseModel=True,
        showMaxMin=True,
        showMode=True,
        showMean=True,
        showMedian=True,
        showEnvelope=True,
        envelopeType=envelope_type,
        showSingle=True,
        singleType=single_type,
        timespan=timespan,
        pcolor=pqlx,
    )
    assert rval is None
    if f_name is None:
        assert os.path.exists("f_name_" + chan + ".png")
        os.remove("f_name_" + chan + ".png")
    else:
        assert os.path.exists(f_name + "_" + chan + ".png")
        os.remove(f_name + "_" + chan + ".png")
    plt.close('all')


@pytest.mark.parametrize("style", ["psd", "pdf"])
@pytest.mark.parametrize("envelope_type", ["10_90", "05_95"])
@pytest.mark.parametrize("f_name", [None, "ftest"])
@pytest.mark.parametrize("single_type", ["5"])
@pytest.mark.parametrize("timespan", [0])
def test_psdPlot_multi(style, envelope_type, f_name, single_type, timespan):
    # data = 'test_data/multi_chan_mseed.mseed'
    data = "test_data/multi_day_multi_channel_AK.mseed"
    st = obspy.read(data)
    rval = psdPlot(
        st,
        style=style,
        f_name=f_name,
        showNoiseModel=True,
        showMaxMin=True,
        showMode=True,
        showMean=True,
        showMedian=True,
        showEnvelope=True,
        envelopeType=envelope_type,
        showSingle=True,
        singleType=single_type,
        timespan=timespan,
        pcolor=pqlx,
    )
    assert rval is None

    sta = str(st[0].stats["station"])
    if f_name is None:
        assert os.path.exists("f_name_" + sta + ".png")
        os.remove("f_name_" + sta + ".png")
    else:
        assert os.path.exists(f_name + "_" + sta + ".png")
        os.remove(f_name + "_" + sta + ".png")
    plt.close('all')


@pytest.mark.parametrize("style", ["psd", "pdf"])
@pytest.mark.parametrize("envelope_type", ["10_90", "05_95"])
@pytest.mark.parametrize("f_name", [None, "ftest"])
@pytest.mark.parametrize("single_type", ["5"])
@pytest.mark.parametrize("timespan", [0])
def test_db_psdPlot_7A(style, envelope_type, f_name, single_type, timespan):
    create_test_db(sql_file="test_data/pycheron_test_db_7A_stream.sql")
    db_path = "test_data/pycheron_test_db.db"
    db = Database(db_path, "WFpych", True)
    chan = "BHE"
    rval = psdPlot(
        db,
        style=style,
        f_name=f_name,
        showNoiseModel=True,
        showMaxMin=True,
        showMode=True,
        showMean=True,
        showMedian=True,
        showEnvelope=True,
        envelopeType=envelope_type,
        showSingle=True,
        singleType=single_type,
        network="7A",
        station="CABN",
        channel=chan,
        timespan=timespan,
        pcolor=pqlx,
    )
    assert rval is None

    if f_name is None:
        assert os.path.exists("f_name_" + chan + ".png")
        os.remove("f_name_" + chan + ".png")
    else:
        assert os.path.exists(f_name + "_" + chan + ".png")
        os.remove(f_name + "_" + chan + ".png")
    plt.close('all')


@pytest.mark.parametrize("style", ["psd", "pdf"])
@pytest.mark.parametrize("envelope_type", ["10_90", "05_95"])
@pytest.mark.parametrize("f_name", [None, "ftest"])
@pytest.mark.parametrize("single_type", ["5"])
@pytest.mark.parametrize("timespan", [0])
def test_db_psdPlot_multi(style, envelope_type, f_name, single_type, timespan):
    create_test_db(sql_file="test_data/multi_day_multi_channel_AK.sql")
    db_path = "test_data/pycheron_test_db.db"
    db = Database(db_path, "WFpych", True)
    sta = "GHO"
    rval = psdPlot(
        db,
        style=style,
        f_name=f_name,
        showNoiseModel=True,
        showMaxMin=True,
        showMode=True,
        showMean=True,
        showMedian=True,
        showEnvelope=True,
        envelopeType=envelope_type,
        showSingle=True,
        singleType=single_type,
        network="AK",
        station=sta,
        timespan=timespan,
        pcolor=pqlx,
    )
    assert rval is None

    if f_name is None:
        assert os.path.exists("f_name_" + sta + ".png")
        os.remove("f_name_" + sta + ".png")
    else:
        assert os.path.exists(f_name + "_" + sta + ".png")
        os.remove(f_name + "_" + sta + ".png")
    plt.close('all')
