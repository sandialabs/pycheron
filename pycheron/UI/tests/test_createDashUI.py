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
import pandas as pd
import pycheron.UI.createDashUI as UI
from pycheron.db.sqllite_db import Database
from pycheron.test.create_test_db import create_test_db
import json

create_test_db(sql_file="test_data/UI-wfdisc-sql.sql")
db_path = "test_data/pycheron_test_db.db"
db = Database(db_path, "UITest", True)


@pytest.fixture(scope="module")
def weight_dict():
    return [
        {
            "min_mask": "",
            "max_mask": "",
            "median_mask": "",
            "variance_mask": "",
            "std_mask": "",
            "rms_mask": "",
            "min_threshold_exceeded": "",
            "max_threshold_exceeded": "",
            "median_threshold_exceeded": "",
            "mean_threshold_exceeded": "",
            "variance_threshold_exceeded": "",
            "std_threshold_exceeded": "",
            "highAmp": "",
            "pegged": "",
            "highAmp_rms_mask": "",
            "pegged_mask": "",
            "num_cals_detected": "",
            "daily_dc_offset_value": "",
            "dc_offset_times": "",
            "dead_chan_mean_metric_masks": "",
            "dead_chan_ADF_metric_masks": "",
            "is_dead_channel": "",
            "total_gaps": "",
            "total_overlaps": "",
            "dc_mask": "",
            "low_amp_mask": "",
            "noise1_mask": "",
            "noise2_mask": "",
            "hi_amp_mask": "",
            "bad_resp_mask": "",
            "dead_chan_exp_hourly_masks": "",
            "dead_chan_lin_hourly_masks": "",
            "dead_chan_gsn_hourly_masks": "",
            "repAmp": "",
            "total_spike_count": "",
            "calibration_signal_counts": "",
            "event_begin_counts": "",
            "event_end_counts": "",
            "event_in_progress_counts": "",
            "negative_leap_counts": "",
            "positive_leap_counts": "",
            "time_correction_applied_counts": "",
            "amplifier_saturation_counts": "",
            "digital_filter_charging_counts": "",
            "digitizer_clipping_counts": "",
            "glitches_counts": "",
            "missing_padded_data_counts": "",
            "spikes_counts": "",
            "suspect_time_tag_counts": "",
            "telemetry_sync_error_counts": "",
            "clock_locked_counts": "",
            "end_time_series_counts": "",
            "long_record_read_counts": "",
            "short_record_read_counts": "",
            "start_time_series_counts": "",
            "station_volume_counts": "",
            "record_count": "",
            "num_records_used": "",
            "noTime": "",
            "poorTQ": "",
            "suspectTime": "",
            "ampSat": "",
            "digFilterChg": "",
            "clip": "",
            "spikes": "",
            "glitch": "",
            "missingPad": "",
            "tsyncErrors": "",
            "calib": "",
            "timingCor": "",
            "noTimeMasks": "",
            "poorTQMasks": "",
            "suspectTimeMasks": "",
            "ampSatMasks": "",
            "digFilterChgMasks": "",
            "clipMasks": "",
            "spikesMasks": "",
            "glitchMasks": "",
            "missingPadMasks": "",
            "tsyncMasks": "",
            "calibMasks": "",
            "tcMasks": "",
            "noise_masks": "",
            "microseism_masks": "",
            "banded_masks": "",
            "snr_masks": "",
        }
    ]


@pytest.mark.parametrize("tab", ["tab-1", "tab-2", "tab-3", "tab-4", "not-tab"])
def test_createDashUI(tab):
    generated_tab = UI._render_content(tab)
    assert generated_tab and isinstance(generated_tab, str)


@pytest.mark.parametrize("clicks", [0, 1, 5])
def test_connect(clicks):
    num_clicks = UI._connect(clicks)
    assert num_clicks is not None and str(clicks) in num_clicks


def test_store_data():
    clicks = str("""{"n_clicks": 1}""")
    map_data = UI._store_data(clicks, "test_data/pycheron_test_db.db")
    assert map_data is not None and isinstance(map_data, str)


def test_gen_map():
    map_data = UI._make_data(db)
    gen_map = UI._gen_map(map_data)
    assert gen_map is not None and isinstance(gen_map, dict)


def test_color_scale():
    map_data = UI._make_data(db)
    color = UI._color_scale(map_data, [0])
    assert color is not None and isinstance(color, list)


@pytest.mark.parametrize("clicks", ["""{"n_clicks": 0}""", """{"n_clicks": 1}"""])
def test_get_network_selections(clicks):
    net_select = UI._get_network_selections(clicks, "test_data/pycheron_test_db.db")
    con = json.loads(net_select)
    if "0" in clicks:
        assert con["response"]["props"]["options"] == []
    else:
        assert con["response"]["props"]["options"] != [] and isinstance(con["response"]["props"]["options"], list)


@pytest.mark.parametrize("clicks", ["""{"n_clicks": 0}""", """{"n_clicks": 1}""", None])
def test_make_map(clicks):
    network = """{"network": "UU"}"""
    map_data = """{"network":{"0":"BW","1":"UU","2":"UU"},
                   "station":{"0":"BGLD","1":"BGU","2":"ZNPU"},
                   "channel":{"0":"EHE","1":"HHE","2":"HHN"},
                   "location":{"0":"--","1":"--","2":"--"},
                   "Overall Quality":{"0":"Bad","1":"Good","2":"Bad"},
                   "latitude":{"0":40.92083,"1":40.92083,"2":40.92083},
                   "longitude":{"0":-113.02983,"1":-113.02983,"2":-113.02983}}"""
    value = UI._make_map(map_data, clicks, network)
    assert value is not None and isinstance(value, str)


@pytest.mark.parametrize("clicks", ["""{"n_clicks": 0}""", """{"n_clicks": 1}""", None])
def test_make_table(clicks):
    network = """{"network": "UU"}"""
    map_data = """{"network":{"0":"BW","1":"UU","2":"UU"},
                   "station":{"0":"BGLD","1":"BGU","2":"ZNPU"},
                   "channel":{"0":"EHE","1":"HHE","2":"HHN"},
                   "location":{"0":"--","1":"--","2":"--"},
                   "Overall Quality":{"0":"Bad","1":"Good","2":"Bad"},
                   "latitude":{"0":40.92083,"1":40.92083,"2":40.92083},
                   "longitude":{"0":-113.02983,"1":-113.02983,"2":-113.02983}}"""
    table = UI._make_table(map_data, clicks, network)
    table = json.loads(table)
    val = table["response"]["props"]["data"]

    if clicks is not None and "1" in clicks:
        assert val is not None and isinstance(val, list)
    else:
        assert val == [{}]


@pytest.mark.parametrize("value", ["Station Ranking", "Other Value"])
def test_disable(value):
    val = UI._disable(value)
    val = json.loads(val)
    chan = val["response"]["network-channel-selector"]["disabled"]
    rank = val["response"]["network-rank-selector"]["disabled"]
    if value == "Station Ranking":
        assert not chan and not rank
    else:
        assert chan and rank


@pytest.mark.parametrize("clicks", ["""{"n_clicks": 0}""", """{"n_clicks": 1}""", None])
def test_get_plot_options(clicks):
    map_data = """{"network":{"0":"BW","1":"UU","2":"UU"},
                   "station":{"0":"BGLD","1":"BGU","2":"ZNPU"},
                   "channel":{"0":"EHE","1":"HHE","2":"HHN"},
                   "location":{"0":"--","1":"--","2":"--"},
                   "Overall Quality":{"0":"Bad","1":"Good","2":"Bad"},
                   "latitude":{"0":40.92083,"1":40.92083,"2":40.92083},
                   "longitude":{"0":-113.02983,"1":-113.02983,"2":-113.02983}}"""
    db_path = "test_data/pycheron_test_db.db"
    plot_opts = UI._get_plot_options(map_data, map_data, clicks, db_path)
    assert plot_opts is not None


@pytest.mark.parametrize("clicks", ["""{"n_clicks": 0}""", """{"n_clicks": 1}""", None])
def test_get_stations_selections(clicks):
    network = """{"network":{"0":"BW","1":"UU","2":"UU"}}"""
    value = "test_data/pycheron_test_db.db"
    sta_sel = UI._get_station_selections(network, clicks, value)
    sta_sel = json.loads(sta_sel)
    assert sta_sel is not None and isinstance(sta_sel, dict)


@pytest.mark.parametrize("clicks, metric", [("""{"n_clicks": 1}""", "deadChannelMetric"), (None, None)])
def test_top_metric_table_plot(clicks, metric):
    network = """{"network":"UU"}"""
    station = """{"station":"DUG"}"""
    channel = """{"channel":"EHE"}"""
    path = "test_data/pycheron_test_db.db"
    data = UI._top_metric_table_plot(network, station, channel, clicks, metric, path)
    assert data is not None


@pytest.mark.parametrize("clicks", ["""{"n_clicks": 0}""", """{"n_clicks": 1}""", None])
def test_get_channel_selections(clicks):
    path = "test_data/pycheron_test_db.db"
    network = """{"network":"UU"}"""
    station = """{"station":"DUG"}"""
    chan_select = UI._get_channel_selections(network, station, clicks, path)
    chan_select = json.loads(chan_select)
    val = chan_select["response"]["props"]["options"]
    if clicks is not None and "1" in clicks:
        assert isinstance(val, list) and val != []
    else:
        assert val == []


@pytest.mark.parametrize("clicks", ["""{"n_clicks": 0}""", """{"n_clicks": 1}""", None])
def test_get_metric_top_panel_selections(clicks):
    network = """{"network":"UU"}"""
    station = """{"station":"DUG"}"""
    channel = """{"channel":"EHE"}"""
    path = "test_data/pycheron_test_db.db"
    top_panel = UI._get_metric_top_panel_selections(
        network, station, channel, clicks, path
    )
    top_panel = json.loads(top_panel)
    val = top_panel["response"]["props"]["options"]
    if clicks is not None and "1" in clicks:
        assert isinstance(val, list) and val[0]["value"] == "DCOffsetTimesMetric"
    else:
        assert isinstance(val, list) and val[0]["value"] == ""


@pytest.mark.parametrize("open_clicks", ["""{"n_clicks_timestamp": 0}""", """{"n_clicks_timestamp": 1}"""])
@pytest.mark.parametrize("close_clicks", ["""{"n_clicks_timestamp": 0}""", """{"n_clicks_timestamp": 1}"""])
@pytest.mark.parametrize("csv_clicks", ["""{"n_clicks_timestamp": 0}""", """{"n_clicks_timestamp": 1}"""])
def test_control_modal(open_clicks, close_clicks, csv_clicks, weight_dict):
    create_test_db(sql_file="test_data/UI-wfdisc-sql.sql")
    db_path = "test_data/pycheron_test_db.db"

    sum_tab = [{"network": "UU", "station": "DUG", "channel": "EHE", "all_counts": 0}]

    if (open_clicks > close_clicks) and (open_clicks > csv_clicks):
        modal = UI.control_modal(
            open_clicks,
            close_clicks,
            csv_clicks,
            str(0),
            [],
            [],
            [],
            db_path,
            [],
            [],
            weight_dict,
        )
        json_modal = json.loads(modal)
        assert len(json_modal["response"]["summary-report-table-modal"]["data"]) > 0
    if (close_clicks > open_clicks) and (close_clicks > csv_clicks):
        modal = UI.control_modal(
            open_clicks,
            close_clicks,
            csv_clicks,
            str(0),
            [],
            [],
            [],
            db_path,
            [],
            [],
            weight_dict,
        )
        json_modal = json.loads(modal)
        assert len(json_modal["response"]["summary-report-table-modal"]["data"]) == 0
    if (csv_clicks > open_clicks) and (csv_clicks > close_clicks):
        modal = UI.control_modal(
            open_clicks,
            close_clicks,
            csv_clicks,
            str(0),
            [],
            [],
            [],
            db_path,
            sum_tab,
            [],
            weight_dict,
        )
        json_modal = json.loads(modal)
        assert len(json_modal["response"]["summary-report-table-modal"]["data"]) == 0

        assert os.path.exists("output.xlsx")
        os.remove("output.xlsx")


def test_summary_report_query():
    create_test_db(sql_file="test_data/UI-wfdisc-sql.sql")
    db_path = "test_data/pycheron_test_db.db"
    db = Database(db_path, "WFpych", True)
    (
        full_table,
        sums_table,
        values_table,
        weights_frame_display,
    ) = UI._summary_report_query(db)

    assert isinstance(full_table, pd.core.frame.DataFrame)
    assert isinstance(sums_table, pd.core.frame.DataFrame)
    assert isinstance(values_table, pd.core.frame.DataFrame)
    assert isinstance(weights_frame_display, pd.core.frame.DataFrame)

    assert full_table.empty is False
    assert sums_table.empty is False
    assert values_table.empty is False
    assert weights_frame_display.empty is False

    assert len(full_table["all_counts"]) > 0
    assert len(sums_table["all_counts"]) > 0
    assert (full_table["all_counts"] == sums_table["all_counts"]).all()

    assert (full_table["range"] == sums_table["range"]).all()
    assert (full_table["range"] == values_table["range"]).all()

    "all_counts" not in values_table.columns


@pytest.mark.parametrize("count", [0, 1, 2, 3, 4])
def test_quality_color(count):
    qcolor = UI._quality_color(count)
    assert isinstance(qcolor, str) and qcolor in ("Good", "Marginal", "Bad", "--")

    if count == 0 or count == 1:
        assert qcolor == "Good"
    if count == 2:
        assert qcolor == "Marginal"
    if count == 3 or count == 4:
        assert qcolor == "Bad"


@pytest.mark.parametrize("count", [-1, None, "TEST"])
def test_quality_color_error(count):
    with pytest.raises(ValueError):
        UI._quality_color(count)


def test_conditional_quality_color():
    COLOR_SCHEMES = {
        "summary_colors": {
            "high_quality": "rgb(92, 184, 92)",
            "medium_quality": "rgb(240, 173, 78)",
            "low_quality": "rgb(217, 83, 79)",
            "value_shade": "rgb(176, 214, 255)",
        },
    }

    df = pd.DataFrame({"all_counts": [0, 1, 2, 3, 4]})
    qcond_color = UI._conditional_quality_color(df, COLOR_SCHEMES, "test")

    assert len(qcond_color) == len(df["all_counts"])

    assert qcond_color[0]["backgroundColor"] == COLOR_SCHEMES["summary_colors"]["high_quality"]
    assert qcond_color[1]["backgroundColor"] == COLOR_SCHEMES["summary_colors"]["high_quality"]
    assert qcond_color[2]["backgroundColor"] == COLOR_SCHEMES["summary_colors"]["medium_quality"]
    assert qcond_color[3]["backgroundColor"] == COLOR_SCHEMES["summary_colors"]["low_quality"]
    assert qcond_color[4]["backgroundColor"] == COLOR_SCHEMES["summary_colors"]["low_quality"]


def test_organize_nsc_columns():
    nsc_cols = ["network", "station", "channel"]

    create_test_db(sql_file="test_data/UI-wfdisc-sql.sql")
    db_path = "test_data/pycheron_test_db.db"
    db = Database(db_path, "WFpych", True)

    full_table, _, _, _ = UI._summary_report_query(db)
    organized = UI._organize_nsc_cols(nsc_cols, full_table)

    assert organized.columns[0] == nsc_cols[0]
    assert organized.columns[1] == nsc_cols[1]
    assert organized.columns[2] == nsc_cols[2]


def test_organize_nsc_columns_error():
    nsc_cols = ["a", "b", "c"]

    db_path = "test_data/pycheron_test_db.db"
    db = Database(db_path, "WFpych", True)

    full_table, _, _, _ = UI._summary_report_query(db)

    with pytest.raises(KeyError):
        UI._organize_nsc_cols(nsc_cols, full_table)


def test_choose_correct_date():
    sum_tab = [
        {
            "vals": ["a", "b"],
            "range": ["2021-01-01T00:00:00.000000Z", "2021-01-02T00:00:00.000000Z"],
        },
        {
            "vals": ["c", "d"],
            "range": ["2021-01-01T00:00:00.000000Z", "2021-01-02T00:00:00.000000Z"],
        },
        {
            "vals": ["e", "f"],
            "range": ["2021-01-03T00:00:00.000000Z", "2021-01-04T00:00:00.000000Z"],
        },
        {
            "vals": ["g", "h"],
            "range": ["2021-01-03T00:00:00.000000Z", "2021-01-04T00:00:00.000000Z"],
        },
    ]

    df = UI._choose_correct_date(sum_tab, "2021-01-01", "2021-01-02")

    pd.testing.assert_series_equal(df["vals"], pd.Series([["a", "b"], ["c", "d"]]), check_names=False)
