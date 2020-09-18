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

import pytest
import pycheron.UI.createDashUI as UI
from pycheron.db.sqllite_db import Database
import json

db = Database("./test_data/test-UI.db", "UITest", True)


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
    map_data = UI._store_data(clicks, "test_data/test-UI.db")
    assert map_data is not None and isinstance(map_data, str)


def test_gen_map():
    map_data = UI._make_data(db)
    gen_map = UI._gen_map(map_data)
    assert gen_map is not None and isinstance(gen_map, dict)


def test_color_scale():
    map_data = UI._make_data(db)
    color = UI._color_scale(map_data, [1])
    assert color is not None and isinstance(color, list)


@pytest.mark.parametrize("clicks", ["""{"n_clicks": 0}""", """{"n_clicks": 1}"""])
def test_get_network_selections(clicks):
    net_select = UI._get_network_selections(clicks, "test_data/test-UI.db")
    con = json.loads(net_select)
    if "0" in clicks:
        assert con["response"]["props"]["options"] == []
    else:
        assert con["response"]["props"]["options"] != [] and isinstance(
            con["response"]["props"]["options"], list
        )


@pytest.mark.parametrize("clicks", ["""{"n_clicks": 0}""", """{"n_clicks": 1}""", None])
def test_make_map(clicks):
    network = """{"network": "UU"}"""
    map_data = """{"network":{"0":"BW","1":"UU","2":"UU"},
                   "station":{"0":"BGLD","1":"BGU","2":"ZNPU"},
                   "channel":{"0":"EHE","1":"HHE","2":"HHN"},
                   "location":{"0":"--","1":"--","2":"--"},
                   "Overall Quality":{"0":"Low","1":"High","2":"Low"},
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
                   "Overall Quality":{"0":"Low","1":"High","2":"Low"},
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
                   "Overall Quality":{"0":"Low","1":"High","2":"Low"},
                   "latitude":{"0":40.92083,"1":40.92083,"2":40.92083},
                   "longitude":{"0":-113.02983,"1":-113.02983,"2":-113.02983}}"""
    db_path = "test_data/test-UI.db"
    plot_opts = UI._get_plot_options(map_data, map_data, clicks, db_path)
    assert plot_opts is not None


@pytest.mark.parametrize("clicks", ["""{"n_clicks": 0}""", """{"n_clicks": 1}""", None])
def test_get_stations_selections(clicks):
    network = """{"network":{"0":"BW","1":"UU","2":"UU"}}"""
    value = "test_data/test-UI.db"
    sta_sel = UI._get_station_selections(network, clicks, value)
    sta_sel = json.loads(sta_sel)
    assert sta_sel is not None and isinstance(sta_sel, dict)


@pytest.mark.parametrize(
    "clicks, metric", [("""{"n_clicks": 1}""", "deadChannelMetric"), (None, None)]
)
def test_top_metric_table_plot(clicks, metric):
    network = """{"network":"UU"}"""
    station = """{"station":"BGU"}"""
    channel = """{"channel":"HHE"}"""
    path = "test_data/test-UI.db"
    data = UI._top_metric_table_plot(network, station, channel, clicks, metric, path)
    assert data is not None


@pytest.mark.parametrize("clicks", ["""{"n_clicks": 0}""", """{"n_clicks": 1}""", None])
def test_get_channel_selections(clicks):
    path = "test_data/test-UI.db"
    network = """{"network":"UU"}"""
    station = """{"station":"BGU"}"""
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
    station = """{"station":"BGU"}"""
    channel = """{"channel":"HHE"}"""
    path = "test_data/test-UI.db"
    top_panel = UI._get_metric_botttom_panel_selections(
        network, station, channel, clicks, path
    )
    top_panel = json.loads(top_panel)
    val = top_panel["response"]["props"]["options"]
    if clicks is not None and "1" in clicks:
        assert isinstance(val, list) and val[0]["value"] == "basicStatsMetric"
    else:
        assert isinstance(val, list) and val[0]["value"] == ""
