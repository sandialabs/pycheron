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

from pycheron.db.sqllite_db import Database
from pycheron.metrics.tests.utils import get_stream_from_client, get_stream_from_file
import sqlite3
from pycheron.metrics.crossCorrMetric import crossCorrMetric
from pycheron.metrics.correlationMetric import correlationMetric
from pycheron.metrics.gapMetric import gapMetric
from pycheron.metrics.sohMetric import sohMetric
from pycheron.psd.noise.networkNoiseModel import networkNoiseModel
from pycheron.psd.noise.stationNoiseModel import stationNoiseModel
from pycheron.metrics.basicStatsMetric import basicStatsMetric
from pycheron.metrics.calibration import calibrationMetric
from pycheron.metrics.psdMetric import psdMetric
from pycheron.plotting.psdPlot import psdPlot
import os

"""
Issues:
1) Inserting metric with override fails in certain circumstances
2) Fortran code in psd module is unreliable and will seg fault intermentiantly
3) Given a stream to psd plot inserting into a Database fails, when trying to find the metric name
4) to_dataFrame returns an error
"""


db = Database("test.db", overwrite=True)
db1 = Database("gapMetric.db", session_name="test", overwrite=True)
db_psd = Database("psdMetric.db", session_name="test", overwrite=True)


def test_sqlite_db_instantiate():
    assert db is not None


def test_sql_db_getters_and_setters():
    db.set_overwrite(True)
    assert db._overwrite

    db.set_session_name("test")
    assert db._session_name == "test"

    con = db.get_connection()
    assert isinstance(con, sqlite3.Connection)


def calculate_cross_corr_metric():
    st1 = get_stream_from_client(
        "NM", "SLM", "00", "BHZ", "2013-11-12T07:09:45.000", 600
    )
    st2 = get_stream_from_client("NM", "SLM", "00", "BHZ", "2013-11-12T07:09:48", 600)
    tr1 = st1[0]
    tr2 = st2[0]
    return crossCorrMetric(tr1, tr2)


def calculate_correlation_metric():
    st = get_stream_from_client(
        "IU", "ANMO", "00", "BH*", "2013-03-01T00:00:00.000", 1440 * 60
    )
    return correlationMetric(st[0], st[1])


def calculate_gap_metric():
    stream = get_stream_from_file("test_data/sglchan_multgap_multoverlap.mseed")
    return gapMetric(stream)


def calculate_soh_metric():
    st = get_stream_from_file("test_data/qualityflags.mseed")
    return sohMetric(st)


def calculate_network_noise_model():
    stream = get_stream_from_client(
        "AK", "GHO", "", "BNZ", "2012-12-12T00:00:00.000", 60 * 60 * 24 * 3
    )
    return networkNoiseModel(stream)


def calculate_station_noise_model():
    stream = get_stream_from_client(
        "AK", "GHO", "", "BNZ", "2012-12-12T00:00:00.000", 60 * 60 * 24 * 3
    )
    return stationNoiseModel(stream)


def calculate_basic_stats_metric():
    stream = get_stream_from_file("test_data/BGU_HHE_001.mseed")
    return basicStatsMetric(stream)


def calculate_calibration_metric():
    stream = get_stream_from_file("test_data/ZNPU_HHN_006.mseed")
    return calibrationMetric(stream)


def calculate_psd_metric():
    stream = get_stream_from_file("test_data/BGU_HHE_002.mseed")
    return psdMetric(stream)


def test_insert_metric_cross_corr():
    cross_corr = calculate_cross_corr_metric()
    db.insert_metric(cross_corr)
    val = db.get_metric("crossCorrMetric")
    assert val is not None


def test_insert_metric_correlation():
    corr = calculate_correlation_metric()
    db.set_overwrite(True)
    db.insert_metric(corr)
    val = db.get_metric("correlationMetric")
    assert val is not None

    # Re-insert to test override functionality
    db.insert_metric(corr)
    val = db.get_metric("correlationMetric")
    assert val is not None


def test_insert_metric_gap():
    gap_sum, gap_dets = calculate_gap_metric()
    db1.insert_metric((gap_sum, gap_dets))
    val = db1.get_metric("gapMetricChannel")
    assert val is not None

    # Re-insert to test override functionality
    db1.insert_metric((gap_sum, gap_dets))
    val = db1.get_metric("gapMetricChannel")
    assert val is not None


def test_insert_metric_soh():
    db2 = Database("sohMetric.db", session_name="test", overwrite=True)
    soh = calculate_soh_metric()
    db2.insert_metric(soh)
    val = db2.get_metric("sohMetricActivityFlags")
    assert val is not None

    # Re-insert to test override functionality
    db2.insert_metric(soh)
    val = db2.get_metric("sohMetricActivityFlags")
    assert val is not None


def test_insert_network_noise():
    net_noise = calculate_network_noise_model()
    db.insert_metric(net_noise)
    val = db.get_metric("networkNoiseModel")
    assert val is not None

    # Re-insert to test override functionality
    db.insert_metric(net_noise)
    val = db.get_metric("networkNoiseModel")
    assert val is not None


def test_insert_station_noise_model():
    sta_noise = calculate_station_noise_model()
    db.insert_metric(sta_noise)
    val = db.get_metric("stationNoiseModel")
    assert val is not None

    # Re-insert to test override functionality
    db.insert_metric(sta_noise)
    val = db.get_metric("stationNoiseModel")
    assert val is not None


def test_insert_basic_stats():
    basic_stats = calculate_basic_stats_metric()
    db.insert_metric(basic_stats)
    val = db.get_metric("basicStatsMetric")
    assert val is not None

    # Re-insert to test override functionality
    basic_stats = calculate_basic_stats_metric()
    db.insert_metric(basic_stats)
    val = db.get_metric("basicStatsMetric")
    assert val is not None


def test_insert_calibration():
    cal = calculate_calibration_metric()
    db.insert_metric(cal)
    val = db.get_metric("calibrationMetric")
    assert val is not None


def test_insert_psd():
    psd = calculate_psd_metric()
    db_psd.insert_metric(psd)
    val = db_psd.get_metric("psdMetric")
    assert val is not None

    stream = get_stream_from_file("test_data/BGU_HHE_002.mseed")
    psdPlot(stream, database=db_psd)
    val = db_psd.get_metric("psdPlot")
    assert val is not None

    # Re-insert to test override functionality
    psdPlot(stream, database=db_psd)
    val = db_psd.get_metric("psdPlot")
    assert val is not None


def test_connect():
    conn = db.connect("test.db")
    assert conn is not None and isinstance(conn, sqlite3.Connection)


def test_networks():
    networks = db.networks()
    assert networks is not None and isinstance(networks, list)


def test_stations():
    stations = db.stations()
    assert stations != [] and isinstance(stations, list)


def test_channels():
    channels = db.channels()
    assert channels is not None and isinstance(channels, list)


def test_session():
    session = db.session()
    assert session is not None and isinstance(session, list)


def test_metrics():
    metrics = db.metrics()
    assert metrics != [] and isinstance(metrics, list)


def test_to_csv():
    db1.to_csv("test_data")
    files = os.listdir("test_data")
    assert any([".csv" in f for f in files])
