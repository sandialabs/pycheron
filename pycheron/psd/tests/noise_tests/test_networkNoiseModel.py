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

from pycheron.psd.noise.networkNoiseModel import networkNoiseModel
from pycheron.test.create_test_db import create_test_db
from pycheron.db.sqllite_db import Database
import numpy as np
import obspy
import pytest
import os


def test_networkNoiseModel_types():
    data = "test_data/7a_cabn_bhe.884965.tar.mseed"
    st = obspy.read(data)
    out = networkNoiseModel(st, plot=True, fname="nnm_test")
    assert os.path.exists("nnm_test.png")
    os.remove("nnm_test.png")

    assert isinstance(out["psdStats"]["noise_matrix_frequency"][0][0], np.ndarray)
    assert isinstance(out["psdStats"]["mode"][0][0], np.ndarray)
    assert isinstance(out["enModel"][0], list)
    assert isinstance(out["psdStats"]["snclq"][0][0], str)
    assert isinstance(out["network"], str)
    assert isinstance(out["network_noise_image_path"], str)
    assert isinstance(out["metric_name"], str)
    assert isinstance(out["stations_en"][0], str)
    assert isinstance(out["stations_z"], list)
    assert isinstance(out["zModel"], list)


def test_networkNoiseModel_e():
    data = "test_data/7a_cabn_bhe.884965.tar.mseed"
    st = obspy.read(data)
    out = networkNoiseModel(st, plot=True, fname="nnm_test")

    assert os.path.exists("nnm_test.png")
    os.remove("nnm_test.png")

    # Check calculations
    assert out["psdStats"]["noise_matrix_frequency"][0][0][0] == pytest.approx(
        0.005255602595335718, 0.001
    )
    assert out["psdStats"]["noise_matrix_frequency"][0][0][-1] == pytest.approx(
        19.74029856522165, 0.001
    )

    assert out["psdStats"]["mode"][0][0][0] == pytest.approx(-127.0, 0.001)
    assert out["psdStats"]["mode"][0][0][-1] == pytest.approx(-20.0, 0.001)

    assert out["enModel"][0][0] == pytest.approx(190.2731384004353, 0.001)
    assert out["enModel"][0][-1] == pytest.approx(0.050657795103555045, 0.001)
    assert out["enModel"][1][0] == pytest.approx(-127.0, 0.001)
    assert out["enModel"][1][-1] == pytest.approx(-20.0, 0.001)

    # Check strings
    assert out["psdStats"]["snclq"][0][0] == "7A.CABN..BHE.D"
    assert out["network"] == "7A"
    assert out["network_noise_image_path"] == "nnm_test.png"
    assert out["metric_name"] == "networkNoiseModel"
    assert out["stations_en"][0] == "CABN"

    # Check for empty values
    assert out["stations_z"] == []
    assert out["zModel"] == [[], []]


def test_networkNoiseModel_e_with_db_arg():
    data = "test_data/7a_cabn_bhe.884965.tar.mseed"
    st = obspy.read(data)

    create_test_db(sql_file="test_data/pycheron_test_db_7A_stream.sql")
    db_path = "test_data/pycheron_test_db.db"
    db = Database(db_path, "WFpych", True)

    out = networkNoiseModel(st, plot=True, fname="nnm_test", database=db)

    assert os.path.exists("nnm_test.png")
    os.remove("nnm_test.png")

    # Check calculations
    assert eval(out["psdStats"]["noise_matrix_frequency"][0][0])[0] == pytest.approx(
        0.005255602595335718, 0.001
    )
    assert eval(out["psdStats"]["noise_matrix_frequency"][0][0])[-1] == pytest.approx(
        19.74029856522165, 0.001
    )
    assert eval(out["psdStats"]["mode"][0][0])[0] == pytest.approx(-127.0, 0.001)
    assert eval(out["psdStats"]["mode"][0][0])[-1] == pytest.approx(-20.0, 0.001)

    assert out["enModel"][0][0] == pytest.approx(190.2731384004353, 0.001)
    assert out["enModel"][0][-1] == pytest.approx(0.050657795103555045, 0.001)
    assert out["enModel"][1][0] == pytest.approx(-127.0, 0.001)
    assert out["enModel"][1][-1] == pytest.approx(-20.0, 0.001)

    # Check strings
    assert eval(out["psdStats"]["snclq"][0])[0] == "7A.CABN..BHE.D"
    assert out["network"] == "7A"
    assert out["network_noise_image_path"] == "nnm_test.png"
    assert out["metric_name"] == "networkNoiseModel"
    assert out["stations_en"][0] == "CABN"

    # Check for empty values
    assert out["stations_z"] == []
    assert out["zModel"] == [[], []]


def test_networkNoiseModel_z():
    data = "test_data/ZNPU_HHZ_001.mseed"
    st = obspy.read(data)
    out = networkNoiseModel(st, plot=True, fname="nnm_test")

    assert os.path.exists("nnm_test.png")
    os.remove("nnm_test.png")

    # Check calculations
    assert out["psdStats"]["noise_matrix_frequency"][0][0][0] == pytest.approx(
        0.005255602595335718, 0.001
    )
    assert out["psdStats"]["noise_matrix_frequency"][0][0][-1] == pytest.approx(
        46.95060701207919, 0.001
    )

    assert out["psdStats"]["mode"][0][0][0] == pytest.approx(-162.0, 0.001)
    assert out["psdStats"]["mode"][0][0][-1] == pytest.approx(-136.0, 0.001)

    # Check strings
    assert out["psdStats"]["snclq"][0][0] == "UU.ZNPU..HHZ.M"
    assert out["network"] == "UU"
    assert out["network_noise_image_path"] == "nnm_test.png"
    assert out["metric_name"] == "networkNoiseModel"
    assert out["stations_z"][0] == "ZNPU"

    # Check for Z model
    assert out["zModel"][0][0] == pytest.approx(190.2731384004353, 0.001)
    assert out["zModel"][0][-1] == pytest.approx(0.021298979153618305, 0.001)
    assert out["zModel"][1][0] == pytest.approx(-162.0, 0.001)
    assert out["zModel"][1][-1] == pytest.approx(-136.0, 0.001)

    # Check for empty values
    assert out["stations_en"] == []
    assert out["enModel"] == [[], []]


def test_networkNoiseModel_e_db():
    create_test_db(sql_file="test_data/pycheron_test_db_7A_stream.sql")
    db_path = "test_data/pycheron_test_db.db"
    db = Database(db_path, "WFpych", True)

    out = networkNoiseModel(
        db, plot=True, fname="nnm_test", network="7A", station="CABN", channel="BHE"
    )

    assert os.path.exists("nnm_test.png")
    os.remove("nnm_test.png")

    # Check calculations
    assert eval(out["psdStats"]["noise_matrix_frequency"][0][0])[0] == pytest.approx(
        0.005255602595335718, 0.001
    )
    assert eval(out["psdStats"]["noise_matrix_frequency"][0][0])[-1] == pytest.approx(
        19.74029856522165, 0.001
    )
    assert eval(out["psdStats"]["mode"][0][0])[0] == pytest.approx(-127.0, 0.001)
    assert eval(out["psdStats"]["mode"][0][0])[-1] == pytest.approx(-20.0, 0.001)

    assert out["enModel"][0][0] == pytest.approx(190.2731384004353, 0.001)
    assert out["enModel"][0][-1] == pytest.approx(0.050657795103555045, 0.001)
    assert out["enModel"][1][0] == pytest.approx(-127.0, 0.001)
    assert out["enModel"][1][-1] == pytest.approx(-20.0, 0.001)

    # Check strings
    assert eval(out["psdStats"]["snclq"][0])[0] == "7A.CABN..BHE.D"
    assert out["network"] == "7A"
    assert out["network_noise_image_path"] == "nnm_test.png"
    assert out["metric_name"] == "networkNoiseModel"
    assert out["stations_en"][0] == "CABN"

    # Check for empty values
    assert out["stations_z"] == []
    assert out["zModel"] == [[], []]


def test_networkNoiseModel_z_db():
    create_test_db(sql_file="test_data/pycheron_test_db_ZNPU_HHZ_stream.sql")
    db_path = "test_data/pycheron_test_db.db"
    db = Database(db_path, "WFpych", True)

    out = networkNoiseModel(
        db, plot=True, fname="nnm_test", network="UU", station="ZNPU", channel="HHZ"
    )

    assert os.path.exists("nnm_test.png")
    os.remove("nnm_test.png")

    # Check calculations
    assert eval(out["psdStats"]["noise_matrix_frequency"][0][0])[0] == pytest.approx(
        0.005255602595335718, 0.001
    )
    assert eval(out["psdStats"]["noise_matrix_frequency"][0][0])[-1] == pytest.approx(
        46.95060701207919, 0.001
    )

    assert eval(out["psdStats"]["mode"][0][0])[0] == pytest.approx(-162.0, 0.001)
    assert eval(out["psdStats"]["mode"][0][0])[-1] == pytest.approx(-136.0, 0.001)

    # Check strings
    assert eval(out["psdStats"]["snclq"][0])[0] == "UU.ZNPU..HHZ.M"
    assert out["network"] == "UU"
    assert out["network_noise_image_path"] == "nnm_test.png"
    assert out["metric_name"] == "networkNoiseModel"
    assert out["stations_z"][0] == "ZNPU"

    # Check Z Model
    assert out["zModel"][0][0] == pytest.approx(190.2731384004353, 0.001)
    assert out["zModel"][0][-1] == pytest.approx(0.021298979153618305, 0.001)
    assert out["zModel"][1][0] == pytest.approx(-162.0, 0.001)
    assert out["zModel"][1][-1] == pytest.approx(-136.0, 0.001)

    # Check for empty values
    assert out["stations_en"] == []
    assert out["enModel"] == [[], []]
