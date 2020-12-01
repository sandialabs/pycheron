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

import dash
import dash_html_components as html
import dash_core_components as dcc
import plotly.express as px
import dash_table
import pandas as pd
from pycheron.db.sqllite_db import Database
import random
import plotly.graph_objs as go
import numpy as np
import json
import os
import base64
import matplotlib.pyplot as plt
from dotenv import load_dotenv
from dateutil import parser

from pycheron.psd.noise.stationNoiseModel import (
    create_psds_from_metric_table,
    gen_stats_plot_vals,
    simple_or_envelope,
)

from pycheron.psd.noise.networkNoiseModel import (
    unique_stations_psds_stats,
    calc_model_vals,
)

from pycheron.plotting.rankingPlot import (
    calc_stats_from_psds_rankplot,
    calc_power_period_rankplot,
    drop_and_rename_rankplot,
)

from pycheron.plotting.dailyPdfPlot import (
    plot_grid_data_fill_in,
    _find_nearest,
)

from dash.dependencies import Input, Output, State
from plotly.subplots import make_subplots
from pycheron.plotting.basicStatPlots import get_basic_stats_data
from pycheron.plotting.psdPlot import get_psd_plot_data
from pycheron.psd.noise.noiseModel import noiseModel
from pycheron.plotting.dailyPdfPlot import get_pdf_plot_data, plot_grid_data_fill_in

load_dotenv()

load_dotenv()

app = dash.Dash(__name__)

app.config["suppress_callback_exceptions"] = True

mapbox_access_token = os.environ.get("MAPBOX_ACCESS_TOKEN")

COLOR_SCHEMES = {
    "psd": [
        [0.0, "rgb(255,64,226)"],
        [0.125, "rgb(0,0,200)"],
        [0.25, "rgb(0,25,255)"],
        [0.375, "rgb(0,152,255)"],
        [0.5, "rgb(44,255,150)"],
        [0.625, "rgb(151,255,0)"],
        [0.75, "rgb(255,234,0)"],
        [0.875, "rgb(255,111,0)"],
        [1, "rgb(255,0,0)"],
    ],
    "daily_pdf": [
        [0.0, "rgb(0, 0, 255)"],
        [0.1, "rgb(0, 0, 255)"],
        [0.1, "rgb(51, 102, 255)"],
        [0.2, "rgb(51, 102, 255)"],
        [0.2, "rgb(102, 204, 255)"],
        [0.3, "rgb(102, 204, 255)"],
        [0.3, "rgb(102, 255, 204)"],
        [0.4, "rgb(102, 255, 204)"],
        [0.4, "rgb(204, 255, 102)"],
        [0.5, "rgb(204, 255, 102)"],
        [0.5, "rgb(255,165,0)"],
        [0.6, "rgb(255,165,0)"],
        [0.6, "rgb(255, 0, 0)"],
        [0.7, "rgb(255, 0, 0)"],
        [0.7, "rgb(153, 51, 102)"],
        [0.8, "rgb(153, 51, 102)"],
        [0.8, "rgb(102, 9, 32)"],
        [0.9, "rgb(102, 9, 32)"],
        [0.9, "rgb(0, 0, 0)"],
        [1.0, "rgb(0, 0, 0)"],
    ],
}

app.layout = html.Div(
    [
        # Header
        html.Div(
            [
                html.Span(
                    "Pycheron Data Quality Control Summary", className="app-title"
                ),
                html.Span(
                    [
                        dcc.Input(
                            id="input-db",
                            type="text",
                            style={"float": "right", "margin-top": ".5%"},
                        ),
                        html.Button(
                            "Connect",
                            id="button",
                            style={
                                "float": "right",
                                "background-color": "#e7e7e7",
                                "margin-right": "5px",
                                "margin-top": ".5%",
                            },
                            n_clicks=0,
                        ),
                    ]
                ),
                html.Span(
                    "Path to Database:",
                    style={
                        "float": "right",
                        "padding-top": "1.25%",
                        "padding-right": "1%",
                    },
                ),
            ],
            className="row header",
        ),
        # tabs
        dcc.Tabs(
            id="tabs-with-classes",
            value="tab-1",
            parent_className="custom-tabs",
            className="custom-tabs-container",
            children=[
                dcc.Tab(
                    label="Overview",
                    value="tab-1",
                    className="custom-tab",
                    selected_className="custom-tab--selected",
                ),
                dcc.Tab(
                    label="Network",
                    value="tab-2",
                    className="custom-tab",
                    selected_className="custom-tab--selected",
                ),
                dcc.Tab(
                    label="Station",
                    value="tab-3",
                    className="custom-tab",
                    selected_className="custom-tab--selected",
                ),
                dcc.Tab(
                    label="Channel",
                    value="tab-4",
                    className="custom-tab",
                    selected_className="custom-tab--selected",
                ),
            ],
        ),
        html.Div(id="tabs-content-classes"),
        # Storage containers for storing values across tabs
        html.Div(id="clicks", style={"display": "none"}),
        html.Div(id="network", style={"display": "none"}),
        html.Div(id="station", style={"display": "none"}),
        html.Div(id="channel", style={"display": "none"}),
        html.Div(id="map-data", style={"display": "none"}),
    ]
)


# Controls tabs
@app.callback(
    Output("tabs-content-classes", "children"), [Input("tabs-with-classes", "value")]
)
def _render_content(tab):
    if tab == "tab-1":
        return overview_tab
    elif tab == "tab-2":
        return network_tab
    elif tab == "tab-3":
        return station_tab
    elif tab == "tab-4":
        return channel_tab
    else:
        html.Div()


# --------------Overview Tab--------------
overview_tab = html.Div(
    [
        # ---- header -----
        html.H3("Overview"),
        html.Hr(),
        # ---- map -----
        html.Div(
            [
                dcc.Graph(
                    figure=go.Figure(
                        data=[go.Scattermapbox()],
                        layout=go.Layout(
                            font=dict(color="#191A1A"),
                            margin=dict(l=35, r=35, b=35, t=45),
                            hovermode="closest",
                            plot_bgcolor="#fffcfc",
                            legend=dict(font=dict(size=28), orientation="h"),
                            title="Station Map",
                            mapbox=dict(accesstoken=mapbox_access_token, style="light"),
                        ),
                    ),
                    id="map-graph",
                ),
            ],
            className="row",
        ),
        # ---- table -----
        html.Div(
            [
                dash_table.DataTable(
                    id="main-table",
                    data=[
                        {
                            "network": "",
                            "station": "",
                            "channel": "",
                            "location": "",
                            "Overall Quality": "",
                            "latitude": "",
                            "longitude": "",
                        }
                    ],
                    columns=[
                        {"name": i, "id": i}
                        for i in [
                            "network",
                            "station",
                            "channel",
                            "location",
                            "Overall Quality",
                            "latitude",
                            "longitude",
                        ]
                    ],
                    style_table={"maxHeight": "500", "overflowY": "scroll"},
                    n_fixed_rows=1,
                    style_cell_conditional=[
                        {
                            "if": {"row_index": "odd"},
                            "backgroundColor": "rgb(248, 248, 248)",
                        }
                    ],
                    style_cell={
                        "minWidth": "14%",
                        "width": "14%",
                        "maxWidth": "14%",
                        "font_family": "Helvetica",
                        "font_size": "12px",
                        "text-align": "left",
                    },
                    style_header={
                        "font_family": "Helvetica",
                        "font_size": "15px",
                        "text-align": "center",
                    },
                    filtering=True,
                    sorting=True,
                    sorting_type="multi",
                )
            ],
            style={"margin-right": ".5%", "margin-left": "3%"},
        ),
    ],
    className="twelve columns",
)

# --------------Network Tab--------------
network_tab = html.Div(
    [
        html.Div(
            [
                # ---- header -----
                html.Span("Network", style={"font-size": "3.0rem"}),
                html.Span(
                    [
                        dcc.Dropdown(
                            id="network-selector",
                            style={
                                "float": "right",
                                "margin-top": "1%",
                                "width": "300px",
                            },
                            options=[{"label": "", "value": ""}],
                        ),
                    ]
                ),
                html.Span(
                    "Select Network:",
                    style={
                        "float": "right",
                        "margin-top": "1.75%",
                        "padding-right": "10px",
                    },
                ),
                html.Hr(),
            ],
            style={"margin-top": "1%"},
        ),
        # Row containing map and table
        html.Div(
            [  # ---- map -----
                html.Div(
                    [
                        dcc.Graph(
                            figure=go.Figure(
                                data=[go.Scattermapbox()],
                                layout=go.Layout(
                                    font=dict(color="#191A1A"),
                                    margin=dict(l=35, r=35, b=35, t=45),
                                    hovermode="closest",
                                    plot_bgcolor="#fffcfc",
                                    legend=dict(font=dict(size=28), orientation="h"),
                                    title="Station Map",
                                    mapbox=dict(
                                        accesstoken=mapbox_access_token, style="light"
                                    ),
                                ),
                            ),
                            id="map-network",
                        )
                    ],
                    className="six columns",
                ),
                # ---- table -----
                html.Div(
                    [
                        dash_table.DataTable(
                            id="network-table",
                            data=[
                                {
                                    "station": "",
                                    "channel": "",
                                    "location": "",
                                    "Overall Quality": "",
                                    "latitude": "",
                                    "longitude": "",
                                }
                            ],
                            columns=[
                                {"name": i, "id": i}
                                for i in [
                                    "station",
                                    "channel",
                                    "location",
                                    "Overall Quality",
                                    "latitude",
                                    "longitude",
                                ]
                            ],
                            css=[
                                {
                                    "selector": ".dash-cell div.dash-cell-value",
                                    "rule": "display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;",
                                }
                            ],
                            style_table={"maxHeight": "100", "overflowY": "scroll"},
                            n_fixed_rows=1,
                            style_cell_conditional=[
                                {
                                    "if": {"row_index": "odd"},
                                    "backgroundColor": "rgb(248, 248, 248)",
                                }
                            ],
                            style_cell={
                                "minWidth": "10%",
                                "width": "10%",
                                "maxWidth": "10%",
                                "font_family": "Helvetica",
                                "font_size": "12px",
                                "text-align": "left",
                                "textOverflow": "ellipsis",
                            },
                            style_header={
                                "font_family": "Helvetica",
                                "font_size": "15px",
                                "text-align": "center",
                            },
                            filtering=True,
                            sorting=True,
                            sorting_type="multi",
                        )
                    ],
                    className="six columns",
                ),
            ],
            className="row",
        ),
        html.Div(
            [
                # ---- plots ----
                html.H4("Plots"),
                html.Div(
                    [
                        dcc.Dropdown(
                            id="network-plot-selector",
                            style={
                                "float": "left",
                                "margin-top": "0%",
                                "width": "300px",
                                "padding-right": "10px",
                            },
                            options=[{"label": "", "value": ""}],
                            placeholder="Select Plot Type",
                        ),
                        dcc.Dropdown(
                            id="network-channel-selector",
                            style={
                                "float": "left",
                                "margin-top": "0%",
                                "width": "300px",
                                "padding-right": "10px",
                            },
                            options=[{"label": "", "value": ""}],
                            placeholder="Select Channel",
                            disabled=True,
                        ),
                        dcc.Dropdown(
                            id="network-rank-selector",
                            style={
                                "float": "left",
                                "margin-top": "0%",
                                "width": "300px",
                                "padding-right": "10px",
                            },
                            options=[{"label": "", "value": ""}],
                            placeholder="Select Rank Period",
                            disabled=True,
                        ),
                    ]
                ),
                dcc.Loading(
                    id="loading-network-graph",
                    type="default",
                    children=html.Div(
                        [
                            dcc.Graph(id="network-graph", style={"display": "none"}),
                        ],
                    ),
                ),
            ]
        ),
    ],
    className="twelve columns",
)

# --------------Station Tab--------------
station_tab = html.Div(
    [
        html.Div(
            [
                # ---- header -----
                html.Span("Station", style={"font-size": "3.0rem"}),
                html.Span(
                    [
                        dcc.Dropdown(
                            id="station-selector",
                            style={
                                "float": "right",
                                "margin-top": "1%",
                                "width": "300px",
                            },
                            options=[{"label": "", "value": ""}],
                        ),
                    ]
                ),
                html.Span(
                    "Select Station:",
                    style={
                        "float": "right",
                        "margin-top": "1.75%",
                        "padding-right": "10px",
                    },
                ),
                html.Hr(),
            ]
        ),
        html.Div(
            [
                dcc.Dropdown(
                    id="station-plot-selector",
                    style={
                        "float": "left",
                        "margin-top": "0%",
                        "width": "300px",
                        "padding-right": "10px",
                        "margin-bottom": "0%",
                    },
                    options=[{"label": "", "value": ""}],
                    placeholder="Select Plot Type",
                )
            ]
        ),
        html.Div(
            [
                dcc.Loading(
                    id="loading-station-graph",
                    type="default",
                    children=html.Div(
                        [
                            dcc.Graph(
                                id="station-graph",
                                style={"display": "none", "margin-top": "5em"},
                            ),
                        ],
                    ),
                ),
            ],
        ),
        html.Div(
            [
                html.H5(
                    children="", id="station-table-name", style={"margin-top": "7%"}
                ),
                dash_table.DataTable(
                    id="station-table",
                    data=[{"station": "", "channel": "", "location": ""}],
                    columns=[
                        {"name": i, "id": i} for i in ["station", "channel", "location"]
                    ],
                    style_table={
                        "maxHeight": "300",
                        "overflowY": "scroll",
                        "overflowX": "scroll",
                    },
                    style_data={"whiteSpace": "normal"},
                    css=[
                        {
                            "selector": ".dash-cell div.dash-cell-value",
                            "rule": "display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;",
                        }
                    ],
                    n_fixed_rows=1,
                    style_cell_conditional=[
                        {
                            "if": {"row_index": "odd"},
                            "backgroundColor": "rgb(248, 248, 248)",
                        }
                    ],
                    style_cell={
                        "font_family": "Helvetica",
                        "font_size": "12px",
                        "text-align": "left",
                        "minWidth": "150px",
                        "maxWidth": "250px",
                        "whiteSpace": "normal",
                        "textOverflow": "ellipsis",
                    },
                    style_header={
                        "font_family": "Helvetica",
                        "font_size": "15px",
                        "text-align": "center",
                    },
                    filtering=True,
                    sorting=True,
                    sorting_type="multi",
                ),
            ],
            style={"display": "none"},
            id="div-station-table",
            className="row",
        ),
        html.Div(
            [
                html.H5("State of Health Metric Data Quality Flags"),
                dash_table.DataTable(
                    id="soh-table",
                    data=[{"station": "", "channel": "", "location": ""}],
                    columns=[
                        {"name": i, "id": i} for i in ["station", "channel", "location"]
                    ],
                    css=[
                        {
                            "selector": ".dash-cell div.dash-cell-value",
                            "rule": "display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;",
                        }
                    ],
                    style_table={
                        "maxHeight": "500",
                        "overflowY": "scroll",
                        "overflowX": "scroll",
                    },
                    n_fixed_rows=1,
                    style_cell_conditional=[
                        {
                            "if": {"row_index": "odd"},
                            "backgroundColor": "rgb(248, 248, 248)",
                        }
                    ],
                    style_cell={
                        "font_family": "Helvetica",
                        "font_size": "12px",
                        "text-align": "left",
                        "minWidth": "150px",
                        "maxWidth": "250px",
                        "whiteSpace": "normal",
                        "textOverflow": "ellipsis",
                    },
                    style_header={
                        "font_family": "Helvetica",
                        "font_size": "15px",
                        "text-align": "center",
                    },
                    filtering=True,
                    sorting=True,
                    sorting_type="multi",
                ),
            ],
            style={"display": "none"},
            id="div-soh-table",
            className="row",
        ),
        html.Div(
            [
                html.H5("State of Health Metric I/O Clock Flags"),
                dash_table.DataTable(
                    id="soh-table1",
                    data=[{"station": "", "channel": "", "location": ""}],
                    columns=[
                        {"name": i, "id": i} for i in ["station", "channel", "location"]
                    ],
                    css=[
                        {
                            "selector": ".dash-cell div.dash-cell-value",
                            "rule": "display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;",
                        }
                    ],
                    style_table={
                        "maxHeight": "500",
                        "overflowY": "scroll",
                        "overflowX": "scroll",
                    },
                    n_fixed_rows=1,
                    style_cell_conditional=[
                        {
                            "if": {"row_index": "odd"},
                            "backgroundColor": "rgb(248, 248, 248)",
                        }
                    ],
                    style_cell={
                        "font_family": "Helvetica",
                        "font_size": "12px",
                        "text-align": "left",
                        "minWidth": "150px",
                        "maxWidth": "250px",
                        "whiteSpace": "normal",
                        "textOverflow": "ellipsis",
                    },
                    style_header={
                        "font_family": "Helvetica",
                        "font_size": "15px",
                        "text-align": "center",
                    },
                    filtering=True,
                    sorting=True,
                    sorting_type="multi",
                ),
            ],
            style={"display": "none"},
            id="div-soh-table1",
            className="row",
        ),
    ],
    className="twelve columns",
)

# --------------Channel Tab--------------
channel_tab = html.Div(
    [
        # ---- header -----
        html.Span("Channel", style={"font-size": "3.0rem"}),
        html.Span(
            [
                dcc.Dropdown(
                    id="channel-selector",
                    style={"float": "right", "margin-top": "1%", "width": "200px"},
                    options=[{"label": "", "value": ""}],
                ),
            ]
        ),
        html.Span(
            "Select Channel:",
            style={"float": "right", "margin-top": "1.30%", "padding-right": "10px"},
        ),
        html.Hr(),
        # Plot selector - row
        html.Div(
            [
                # ---- Channel Plot Selector Top Panel  -----
                dcc.Dropdown(
                    id="channel-plot-selectors-top-panel",
                    style={
                        "float": "left",
                        "width": "300px",
                        "padding-right": "10px",
                        "margin-bottom": "2%",
                        "margin-top": "0%",
                        "zIndex": "1112",
                        "position": "relative",
                    },
                    options=[{"label": "", "value": ""}],
                    placeholder="Select Plot Type",
                )
            ],
            className="row",
        ),
        # Image container - row, hidden
        html.Div(
            [
                dcc.Loading(
                    id="loading-channel-top-graph",
                    type="default",
                    children=html.Div(
                        [
                            dcc.Graph(id="channel-top-graph"),
                        ],
                    ),
                ),
            ],
            className="row",
            id="div-image-holder-top",
            style={"display": "none"},
        ),
        # Table container - row, hidden
        html.Div(
            [
                dcc.Loading(
                    id="loading-top-metric-table",
                    type="default",
                    children=html.Div(
                        [
                            dash_table.DataTable(
                                id="top-metric-table",
                                data=[
                                    {
                                        "network": "",
                                        "station": "",
                                        "channel": "",
                                        "location": "",
                                        "start_time": "",
                                        "end_time": "",
                                        "dc_offset_times": "",
                                        "masks": "",
                                    }
                                ],
                                columns=[
                                    {"name": i, "id": i}
                                    for i in [
                                        "network",
                                        "station",
                                        "channel",
                                        "location",
                                        "start_time",
                                        "end_time",
                                        "dc_offset_times",
                                        "masks",
                                    ]
                                ],
                                css=[
                                    {
                                        "margin-top": "12%",
                                        "selector": ".dash-cell div.dash-cell-value",
                                        "rule": "display: inline; white-space: inherit; overflow: inherit; "
                                        "text-overflow: inherit;",
                                    }
                                ],
                                n_fixed_rows=1,
                                style_table={
                                    "overflowY": "scroll",
                                    "overflowX": "auto",
                                },
                                style_cell_conditional=[
                                    {
                                        "if": {"row_index": "odd"},
                                        "backgroundColor": "rgb(248, 248, 248)",
                                    }
                                ],
                                style_cell={
                                    "font_family": "Helvetica",
                                    "font_size": "12px",
                                    "text-align": "left",
                                    "minWidth": "150px",
                                    "maxWidth": "250px",
                                    "whiteSpace": "normal",
                                    "textOverflow": "ellipsis",
                                },
                                style_header={
                                    "font_family": "Helvetica",
                                    "font_size": "15px",
                                    "text-align": "center",
                                },
                                filtering=True,
                                sorting=True,
                                sorting_type="multi",
                            ),
                        ],
                    ),
                ),
            ],
            className="row",
            style={"display": "none"},
            id="top-metric-table-div",
        ),
        # ---- Channel Plot Selector Bottom Panel -----
        # Plot selector - row
        html.Div(
            [
                # ---- Channel Plot Selector Bottom Panel  -----
                dcc.Dropdown(
                    id="channel-plot-selectors-bottom-panel",
                    style={
                        "float": "left",
                        "width": "300px",
                        "padding-right": "10px",
                        "margin-bottom": "2%",
                        "margin-top": "1%",
                        "zIndex": "1111",
                        "position": "relative",
                    },
                    options=[{"label": "", "value": ""}],
                    placeholder="Select Plot Type",
                )
            ],
            className="row",
        ),
        # Image container - row, hidden
        html.Div(
            [
                dcc.Loading(
                    id="loading-channel-bottom-graph",
                    type="default",
                    children=html.Div(
                        [
                            dcc.Graph(id="channel-bottom-graph"),
                        ],
                    ),
                ),
            ],
            className="row",
            id="div-image-holder-bottom",
            style={"display": "none", "float": "center", "position": "relative"},
        ),
        # Table container - row, hidden
        html.Div(
            [
                dcc.Loading(
                    id="loading-bottom-metric-table",
                    type="default",
                    children=html.Div(
                        [  # ---
                            dash_table.DataTable(
                                id="bottom-metric-table",
                                data=[
                                    {
                                        "network": "",
                                        "station": "",
                                        "channel": "",
                                        "location": "",
                                        "start_time": "",
                                        "end_time": "",
                                        "dc_offset_times": "",
                                        "masks": "",
                                    }
                                ],
                                columns=[
                                    {"name": i, "id": i}
                                    for i in [
                                        "network",
                                        "station",
                                        "channel",
                                        "location",
                                        "start_time",
                                        "end_time",
                                        "dc_offset_times",
                                        "masks",
                                    ]
                                ],
                                css=[
                                    {
                                        "margin-top": "12%",
                                        "selector": ".dash-cell div.dash-cell-value",
                                        "rule": "display: inline; white-space: inherit; overflow: inherit; "
                                        "text-overflow: inherit;",
                                    }
                                ],
                                n_fixed_rows=1,
                                style_table={
                                    "overflowY": "scroll",
                                    "overflowX": "auto",
                                },
                                style_cell_conditional=[
                                    {
                                        "if": {"row_index": "odd"},
                                        "backgroundColor": "rgb(248, 248, 248)",
                                    }
                                ],
                                style_cell={
                                    "font_family": "Helvetica",
                                    "font_size": "12px",
                                    "text-align": "left",
                                    "minWidth": "150px",
                                    "maxWidth": "250px",
                                    "whiteSpace": "normal",
                                    "textOverflow": "ellipsis",
                                },
                                style_header={
                                    "font_family": "Helvetica",
                                    "font_size": "15px",
                                    "text-align": "center",
                                },
                                filtering=True,
                                sorting=True,
                                sorting_type="multi",
                            ),
                        ]
                    ),
                )
            ],
            className="row",
            style={"display": "none"},
            id="bottom-metric-table-div",
        ),
    ],
    className="twelve columns",
)


# This stores the button clicks across tables
@app.callback(Output("clicks", "children"), [Input("button", "n_clicks")])
def _connect(n_clicks):
    return json.dumps({"n_clicks": n_clicks})


# ------ Map data functions and storage -----------
def _make_data(db):
    pycheron_df = db.view()

    ########## Manipulate data ###########

    # Grab out a subset of columns based on what we need for table/map. Eventually this column list should be
    # 'network', 'station', 'channel', 'location', 'quality', 'lat', 'lon'
    pych_subset_df = pycheron_df.loc[
        :, ["network", "station", "channel", "location", "latitude", "longitude"]
    ]

    # Drop duplicates so just left with unique values
    pych_uniq_df = pych_subset_df.drop_duplicates()

    # Drop rows where value is None, if strings then catch those too
    pych_df = pych_uniq_df.mask(pych_uniq_df.eq("None").dropna())
    pych_Ndf = pych_df.dropna()

    # Drop the correlation metric information from the network column that has 'UU:UU' as the network
    pdf = pych_Ndf[~pych_Ndf["network"].isin(["UU:UU"])]

    # Add in quality column information. Right now this is just a placeholder, but will eventually be part of the database
    # Generate quality list and create emtpy list to append to
    qual_col = []

    for i in range(pdf["network"].count()):
        qual_col.append("--")

    # Make list a dataframe, then append it to pdf pandas dataframe as a column. Make column header before appending.
    qual_pd = pd.DataFrame(qual_col, columns=["Overall Quality"])

    # Sort stations in alphabetical order
    pDF = pdf.sort_values(["station", "channel"], ascending=True)

    # Reset index the rows before appending
    Pdf = pDF.reset_index(drop=True)

    # Join two dataframes together
    map_data = Pdf.join(qual_pd)

    return map_data


@app.callback(
    Output("map-data", "children"),
    [Input("clicks", "children")],
    [State("input-db", "value")],
)
def _store_data(clicks, value):
    clicks_ = json.loads(clicks)
    if clicks_["n_clicks"] > 0:
        db = Database(value)
        map_data = _make_data(db)
        return map_data.to_json()


def _gen_map(map_data):
    # Generate map based on dataframe, use layout_right
    # groupby returns a dictionary mapping the values of the first field
    # 'classification' onto a list of record dictionaries with that
    # classification value.
    map_data["Sta Chan"] = map_data["station"].str.cat(map_data["channel"], sep=" ")
    return {
        "data": [
            {
                "type": "scattermapbox",
                "lat": list(map_data["latitude"]),
                "lon": list(map_data["longitude"]),
                "text": list(map_data["Sta Chan"]),
                "mode": "markers",
                "name": list(map_data["station"]),
                "marker": {"size": 10, "opacity": 1.0, "color": _color_scale(map_data)},
            }
        ],
        "layout": dict(
            autosize=True,
            font=dict(color="#191A1A"),
            titlefont=dict(color="#191A1A", size="18"),
            margin=dict(l=35, r=35, b=35, t=45),
            hovermode="closest",
            plot_bgcolor="#fffcfc",
            legend=dict(font=dict(size=20), orientation="h"),
            title="Station Map",
            mapbox=dict(
                accesstoken=mapbox_access_token,
                style="light",
                center=dict(
                    lat=map_data["latitude"].astype(float).mean(),
                    lon=map_data["longitude"].astype(float).mean(),
                ),
                zoom=3,
            ),
        ),
    }


def _color_scale(md, selected_row_indices=[]):
    color = []
    for row in md["Overall Quality"]:
        # Neon Green
        if row == "High":
            color.append("#26EC04")
        # Neon Yellow
        elif row == "Medium":
            color.append("#F3F315")
        # Neon Red
        elif row == "Low":
            color.append("#ff0101")
        else:
            color.append("#1500FA")
    for i in selected_row_indices:
        color[i] = "#1500FA"
    return color


# ---------- Overview Functions -------------------

# Map
@app.callback(
    Output("map-graph", "figure"),
    [Input("map-data", "children"), Input("clicks", "children")],
)
def _make_map(map_data, clicks):
    clicks_ = json.loads(clicks)
    if clicks_["n_clicks"] > 0:
        df = pd.read_json(map_data)
        return _gen_map(df)
    else:
        return go.Figure(
            data=[go.Scattermapbox()],
            layout=go.Layout(
                font=dict(color="#191A1A"),
                margin=dict(l=35, r=35, b=35, t=45),
                hovermode="closest",
                plot_bgcolor="#fffcfc",
                legend=dict(font=dict(size=28), orientation="h"),
                title="Station Map",
                mapbox=dict(accesstoken=mapbox_access_token, style="light"),
            ),
        )


# Table
@app.callback(
    Output("main-table", "data"),
    [Input("map-data", "children"), Input("clicks", "children")],
)
def _make_table(map_data, clicks):
    clicks_ = json.loads(clicks)
    if clicks_["n_clicks"] > 0:
        df = pd.read_json(map_data)
        return df.to_dict("records")
    else:
        return [{}]


# ---------- Network Functions -------------------

# This gets the networks available
@app.callback(
    Output("network-selector", "options"),
    [Input("clicks", "children")],
    [State("input-db", "value")],
)
def _get_network_selections(children, value):
    clicks = json.loads(children)
    if clicks["n_clicks"] > 0:
        db = Database(value).networks()
        return [{"label": val, "value": val} for val in np.unique(db) if ":" not in val]
    else:
        return []


# This stores the network across tables
@app.callback(Output("network", "children"), [Input("network-selector", "value")])
def _get_network(value):
    return json.dumps({"network": value})


# Map
@app.callback(
    Output("map-network", "figure"),
    [
        Input("map-data", "children"),
        Input("clicks", "children"),
        Input("network", "children"),
    ],
)
def _make_map(map_data, clicks, network):
    if clicks is not None and network is not None:
        clicks_ = json.loads(clicks)
        value = json.loads(network)
        if clicks_["n_clicks"] > 0 and value["network"] is not None:
            load = pd.read_json(map_data)
            df = load.loc[load["network"] == value["network"]]
            return _gen_map(df)
        else:
            return go.Figure(
                data=[go.Scattermapbox()],
                layout=go.Layout(
                    font=dict(color="#191A1A"),
                    margin=dict(l=35, r=35, b=35, t=45),
                    hovermode="closest",
                    plot_bgcolor="#fffcfc",
                    legend=dict(font=dict(size=28), orientation="h"),
                    title="Station Map",
                    mapbox=dict(accesstoken=mapbox_access_token, style="light"),
                ),
            )
    else:
        return go.Figure(
            data=[go.Scattermapbox()],
            layout=go.Layout(
                font=dict(color="#191A1A"),
                margin=dict(l=35, r=35, b=35, t=45),
                hovermode="closest",
                plot_bgcolor="#fffcfc",
                legend=dict(font=dict(size=28), orientation="h"),
                title="Station Map",
                mapbox=dict(accesstoken=mapbox_access_token, style="light"),
            ),
        )


# Table
@app.callback(
    Output("network-table", "data"),
    [
        Input("map-data", "children"),
        Input("clicks", "children"),
        Input("network", "children"),
    ],
)
def _make_table(map_data, clicks, network):
    if clicks is not None and network is not None:
        clicks_ = json.loads(clicks)
        value = json.loads(network)
        if clicks_["n_clicks"] > 0:
            load = pd.read_json(map_data)
            df = load.loc[load["network"] == value["network"]]
            return df.to_dict("records")
        else:
            return [{}]
    else:
        return [{}]


# Plot Selectors
@app.callback(
    [
        Output("network-plot-selector", "options"),
        Output("network-channel-selector", "options"),
        Output("network-rank-selector", "options"),
    ],
    [Input("clicks", "children"), Input("network", "children")],
    [State("input-db", "value")],
)
def _get_plot_options(clicks, network, value):
    if clicks is not None and network is not None:
        clicks_ = json.loads(clicks)
        network_ = json.loads(network)
        connect_ = value
        if clicks_["n_clicks"] > 0:
            # default values for network-channel/rank-selector disabled
            path = connect_.split("/")
            p = ""
            for i in range(1, len(path) - 1):
                p = p + "/" + path[i]
            p = p + "/" + network_["network"]
            network_opt = []
            network_opt.append({"label": "Station Ranking", "value": "Station Ranking"})
            ranks = [1, 6.5, 30, 100]
            rank_opt = [{"label": val, "value": val} for val in ranks]

            db = Database(value)
            unique_nc = (
                db.get_metric("psdMetric")
                .groupby(["network", "channel"])
                .size()
                .reset_index(name="Freq")
            )
            unique_nc = unique_nc.drop("Freq", 1)
            unique_nc = unique_nc.drop("network", 1)
            unique_nc = unique_nc[~unique_nc.isin(["None"]).any(axis=1)]
            unique_collect = unique_nc["channel"].tolist()
            unique_nc_opt = [{"label": val, "value": val} for val in unique_collect]
            stas = ["stas_opt"]
            stas_opt = [{"label": val, "value": val} for val in stas]

            # adding network noise model option
            network_opt.append(
                {"label": "Network Noise Model", "value": "Network Noise Model"}
            )

            return network_opt, unique_nc_opt, rank_opt


# Disablers/Ablers
@app.callback(
    [
        Output("network-channel-selector", "disabled"),
        Output("network-rank-selector", "disabled"),
    ],
    [Input("network-plot-selector", "value")],
)
def _disable(value):
    chan = True
    rank = True
    if value == "Station Ranking":
        chan = False
        rank = False
    return chan, rank


# Plots
@app.callback(
    [
        Output("network-graph", "figure"),
        Output("network-graph", "style"),
    ],
    [
        Input("network-plot-selector", "value"),
        Input("network-channel-selector", "value"),
        Input("network-rank-selector", "value"),
        Input("clicks", "children"),
        Input("network", "children"),
    ],
    [State("input-db", "value")],
)
def _get_plot(plot, chan, rank, clicks, network, path):
    fig = go.Figure()
    style_graph = {"display": "none"}
    if clicks is not None and network is not None and plot is not None:
        clicks_ = json.loads(clicks)
        network_ = json.loads(network)
        connect_ = path
        if clicks_["n_clicks"] > 0:
            # default values for network-channel/rank-selector disabled
            path = connect_.split("/")
            p = ""
            img = ""
            for i in range(1, len(path) - 1):
                p = p + "/" + path[i]
            p = p + "/" + network_["network"]
            db = Database("/".join(path))
            if plot == "Station Ranking" and None not in (chan, rank):
                style_graph = {"display": "block", "margin-top": "5em"}
                psd_channels, stats, period = calc_stats_from_psds_rankplot(
                    database=db, network=network_["network"], channel=chan
                )
                time_start, time_end = (
                    psd_channels[0][0][0][0][-1][-2],
                    psd_channels[0][0][0][0][-1][-1],
                )
                df_dict = calc_power_period_rankplot(stats, period)
                chan_z, chan_e, chan_n = (
                    df_dict["chan_z"],
                    df_dict["chan_e"],
                    df_dict["chan_n"],
                )
                df_z, df_e, df_n = drop_and_rename_rankplot(df_dict)
                if df_z.empty is False:
                    index_dz = _find_nearest(
                        np.asarray(df_z.columns, dtype=float), rank
                    )
                    df_z_rank = df_z.sort_values(by=[df_z.columns[index_dz]])
                    plot_dat = plot_grid_data_fill_in(df_z_rank)
                    ranker = float(df_z.columns[index_dz])
                if df_e.empty is False:
                    index_de = _find_nearest(
                        np.asarray(df_e.columns, dtype=float), rank
                    )
                    df_e_rank = df_e.sort_values(by=[df_e.columns[index_de]])
                    plot_dat = plot_grid_data_fill_in(df_e_rank)
                    ranker = float(df_e.columns[index_de])
                if df_n.empty is False:
                    index_dn = _find_nearest(
                        np.asarray(df_n.columns, dtype=float), rank
                    )
                    df_n_rank = df_n.sort_values(by=[df_n.columns[index_dn]])
                    plot_dat = plot_grid_data_fill_in(df_n_rank)
                    ranker = float(df_n.columns[index_dn])
                plot_dat = plot_dat.reindex(index=plot_dat.index[::-1])
                stations = plot_dat.axes[0].tolist()
                seconds = plot_dat.axes[1].tolist()
                x = [float(sec) for sec in seconds]
                y = stations
                z_vals = []
                for i in range(plot_dat.shape[0]):
                    col_vals = plot_dat.iloc[i]
                    z = []
                    for col_val in col_vals:
                        z.append(col_val)
                    z_vals.append(z)

                height_scale = 300 if len(y) == 1 else 600
                colorbar_scale = 2 if len(y) == 1 else 1
                rank_title = (
                    "Station Ranking Plot for Channel: "
                    + chan
                    + ", Ranked by: "
                    + str(ranker)
                    + "s"
                    + " ("
                    + str(time_start)
                    + " -- "
                    + str(time_end)
                    + ")"
                )
                colormap = COLOR_SCHEMES["daily_pdf"]

                fig = go.Figure(
                    data=go.Heatmap(
                        x=x,
                        y=y,
                        z=z_vals,
                        zmin=-10,
                        zmax=80,
                        zauto=False,
                        colorscale=colormap,
                        colorbar=dict(
                            len=colorbar_scale,
                            title="D",
                            titleside="top",
                            tickmode="array",
                            tickvals=np.linspace(-5, 75, 10),
                            ticktext=[
                                "D<0",
                                "10>D>=0",
                                "20>d>=10",
                                "30>D>=20",
                                "40>D>=30",
                                "50>D>=40",
                                "60>D>=50",
                                "70>D>=60",
                                "D>=70",
                                "No Data",
                            ],
                            ticks="outside",
                        ),
                    )
                )
                # Correct size of bounding box for log scale
                fig.add_shape(
                    type="rect",
                    xref="x domain",
                    yref="y domain",
                    x0=ranker - ranker / 24,
                    x1=ranker + ranker / 24,
                    y0=-0.5,
                    y1=len(y) - 0.5,
                    line=dict(color="Black"),
                )
                fig.update_xaxes(type="log")

                axis_template_y = dict(
                    showgrid=True,
                    zeroline=True,
                    linecolor="black",
                    showticklabels=True,
                    ticks="outside",
                    dtick=1,
                )
                axis_template_x = dict(
                    showgrid=True,
                    zeroline=True,
                    linecolor="black",
                    showticklabels=True,
                    ticks="outside",
                    tickvals=x,
                )

                fig.update_layout(
                    xaxis=axis_template_x,
                    yaxis=axis_template_y,
                    xaxis_title="Period (seconds)",
                    yaxis_title="Station",
                    title=rank_title,
                    height=height_scale,
                    width=1600,
                )

            if plot == "Network Noise Model":
                net_title = "Noise Model for Network: " + network_["network"]
                style_graph = {"display": "block", "margin-top": "5em"}
                df_dict = unique_stations_psds_stats(db, network_["network"])
                out = calc_model_vals(df_dict, network_["network"])
                if len(out["enModel"][0]) > 0 and len(out["zModel"][0]) == 0:
                    d = {
                        "Period": np.log10(out["enModel"][0]),
                        "Power (dB)": out["enModel"][1],
                        "label": ["E/N Network Model"] * len(out["enModel"][0]),
                    }
                    df = pd.DataFrame(data=d)
                elif len(out["zModel"][0]) > 0 and len(out["enModel"][0]) == 0:
                    d = {
                        "Period": np.log10(out["zModel"][0]),
                        "Power (dB)": out["zModel"][1],
                        "label": ["Z Network Model"] * len(out["zModel"][0]),
                    }
                    df = pd.DataFrame(data=d)
                else:
                    den = {
                        "Period": np.log10(out["enModel"][0]),
                        "Power (dB)": out["enModel"][1],
                        "label": ["E/N Network Model"] * len(out["enModel"][0]),
                    }
                    dfen = pd.DataFrame(data=den)
                    dz = {
                        "Period": np.log10(out["zModel"][0]),
                        "Power (dB)": out["zModel"][1],
                        "label": ["Z Network Model"] * len(out["zModel"][0]),
                    }
                    dfz = pd.DataFrame(data=dz)
                    df = pd.concat([dfen, dfz])
                fig = px.line(
                    df,
                    x="Period",
                    y="Power (dB)",
                    color="label",
                    labels={"Period": "Period (10^n)"},
                    title=net_title,
                    height=600,
                    width=1200,
                )
    return (fig, style_graph)


# ---------- Station Functions -------------------

# This gets the station available for specific network
@app.callback(
    Output("station-selector", "options"),
    [Input("network", "children"), Input("clicks", "children")],
    [State("input-db", "value")],
)
def _get_station_selections(network, clicks, value):
    if clicks is not None and network is not None:
        clicks_ = json.loads(clicks)
        network_ = json.loads(network)
        if clicks_["n_clicks"] > 0:
            db = Database(value).view(network=network_["network"]).station
            return [
                {"label": val, "value": val} for val in np.unique(db) if ":" not in val
            ]
        else:
            return []
    else:
        return []


# This stores the station across tables
@app.callback(Output("station", "children"), [Input("station-selector", "value")])
def _get_station(value):
    return json.dumps({"station": value})


@app.callback(
    Output("station-plot-selector", "options"),
    [
        Input("network", "children"),
        Input("station", "children"),
        Input("clicks", "children"),
    ],
    [State("input-db", "value")],
)
def _get_plot_options(network, station, clicks, value):
    if clicks is not None and network is not None and station is not None:
        clicks_ = json.loads(clicks)
        network_ = json.loads(network)
        if clicks_["n_clicks"] > 0:
            station_ = json.loads(station)
            db = Database(value)
            opt = list(
                np.unique(
                    db.view(
                        network=network_["network"],
                        station=station_["station"],
                        channel="None",
                    ).metric
                )
            )
            noise = list(
                np.unique(
                    db.view(
                        metric_name="stationNoiseModel", station=station_["station"]
                    ).metric
                )
            )
            if noise:
                opt.append(noise[0])

            return [{"label": val, "value": val} for val in opt]
        else:
            return [{"label": "", "value": ""}]
    else:
        return [{"label": "", "value": ""}]


@app.callback(
    [
        Output("station-plot", "style"),
        Output("div-station-table", "style"),
        Output("div-soh-table", "style"),
        Output("div-soh-table1", "style"),
        Output("station-plot", "src"),
        Output("station-table", "data"),
        Output("station-table", "columns"),
        Output("soh-table", "data"),
        Output("soh-table", "columns"),
        Output("soh-table1", "data"),
        Output("soh-table1", "columns"),
        Output("station-table-name", "children"),
        Output("station-graph", "figure"),
        Output("station-graph", "style"),
    ],
    [
        Input("network", "children"),
        Input("clicks", "children"),
        Input("station", "children"),
        Input("station-plot-selector", "value"),
    ],
    [State("input-db", "value")],
)
def _get_plot(network, clicks, station, plot, db_path):
    if clicks is not None and network is not None and station is not None:
        clicks_ = json.loads(clicks)
        network_ = json.loads(network)
        connect_ = db_path
        # Initializing values
        style_img = {"display": "none"}
        style_table = {"display": "none"}
        style_soh = {"display": "none"}
        # For gapMetricStation and sohMetricActivity
        tb = [{"station": "", "channel": "", "location": ""}]
        cols = [{"name": i, "id": i} for i in ["station", "channel", "location"]]
        # For sohMetricDataQualityFlags
        tb1 = [{"station": "", "channel": "", "location": ""}]
        cols1 = [{"name": i, "id": i} for i in ["station", "channel", "location"]]
        # For sohMetricIOClockFlags
        tb2 = [{"station": "", "channel": "", "location": ""}]
        cols2 = [{"name": i, "id": i} for i in ["station", "channel", "location"]]
        src = ""
        name = ""
        style_graph = {"display": "none", "margin-top": "5em"}
        if clicks_["n_clicks"] > 0:
            station_ = json.loads(station)

            path = connect_.split("/")
            p = ""
            for i in range(1, len(path) - 1):
                p = p + "/" + path[i]
            p = p + "/" + network_["network"] + "/" + station_["station"]
            default_session = "WFpych"
            db = Database(db_path)
            fig = go.Figure()

            if plot == "stationNoiseModel":
                style_graph = {"display": "block", "margin-top": "5em"}
                db_metrics = db.get_metric(
                    "stationNoiseModel",
                    network_["network"],
                    station_["station"],
                    session=default_session,
                )
                # TODO: Check how to grab specific type/location. Otherwise, go with default
                psd_type = db_metrics["type"].iloc[0]
                sta_title = (
                    "Noise Model for Station: "
                    + station_["station"]
                    + " ("
                    + psd_type
                    + ")"
                )

                psds = create_psds_from_metric_table(
                    db,
                    network_["network"],
                    station_["station"],
                    session=default_session,
                )
                df_dict = gen_stats_plot_vals(psds, type=psd_type)
                out = simple_or_envelope(
                    df_dict, network_["network"], station_["station"], type=psd_type
                )
                if psd_type.startswith("s"):
                    e_sta = out["e_station_noiseModel"].to_frame()
                    e_sta["sta"] = ["E station"] * e_sta.shape[0]

                    n_sta = out["n_station_noiseModel"].to_frame()
                    n_sta["sta"] = ["N station"] * n_sta.shape[0]

                    z_sta = out["z_station_noiseModel"].to_frame()
                    z_sta["sta"] = ["Z station"] * z_sta.shape[0]

                    final_df = pd.concat([e_sta, n_sta, z_sta])
                    final_df["Period"] = np.log10(final_df.index)
                    final_df.columns = ["Power(dB)", "sta", "Period"]
                    fig = px.line(final_df, x="period", y="Power(dB)", color="sta")
                else:
                    e_low = out["e_station_noiseModel_low"].to_frame()
                    e_low["sta"] = ["E station"] * e_low.shape[0]
                    e_low["hl"] = ["low"] * e_low.shape[0]

                    e_high = out["e_station_noiseModel_high"].to_frame()
                    e_high["sta"] = ["E station"] * e_high.shape[0]
                    e_high["hl"] = ["high"] * e_high.shape[0]

                    n_low = out["n_station_noiseModel_low"].to_frame()
                    n_low["sta"] = ["N station"] * n_low.shape[0]
                    n_low["hl"] = ["low"] * n_low.shape[0]

                    n_high = out["n_station_noiseModel_high"].to_frame()
                    n_high["sta"] = ["N station"] * n_high.shape[0]
                    n_high["hl"] = ["high"] * n_high.shape[0]

                    z_low = out["z_station_noiseModel_low"].to_frame()
                    z_low["sta"] = ["Z station"] * z_low.shape[0]
                    z_low["hl"] = ["low"] * z_low.shape[0]

                    z_high = out["z_station_noiseModel_high"].to_frame()
                    z_high["sta"] = ["Z station"] * z_high.shape[0]
                    z_high["hl"] = ["high"] * z_high.shape[0]

                    final_df = pd.concat([e_low, e_high, n_low, n_high, z_low, z_high])
                    final_df["Period"] = np.log10(final_df.index)
                    final_df.columns = ["Power(dB)", "sta", "hl", "Period"]
                    fig = px.line(
                        final_df,
                        x="Period",
                        y="Power(dB)",
                        color="sta",
                        line_group="hl",
                        hover_name="hl",
                        labels={"Period": "Period (10^n)"},
                        title=sta_title,
                        height=600,
                        width=1200,
                    )

            # return table
            if plot == "gapMetricStation" or plot == "sohMetric":
                style_table = {"display": "block"}
                db = Database(db_path)
                if plot == "sohMetric":
                    name = "State of Health Metric Activity Flags"
                    style_soh = {"display": "block"}
                    # Table 1, sohActivity
                    af = db.get_metric(
                        metric_name=plot + "ActivityFlags",
                        network=network_["network"],
                        station=station_["station"],
                    )
                    tb = af.to_dict("records")
                    cols = [{"name": i, "id": i} for i in af.columns]
                    # Table 2, sohDQ
                    df = db.get_metric(
                        metric_name=plot + "DataQualityFlags",
                        network=network_["network"],
                        station=station_["station"],
                    )
                    tb1 = df.to_dict("records")
                    cols1 = [{"name": i, "id": i} for i in df.columns]

                    # Table 3, IOClock
                    iof = db.get_metric(
                        metric_name=plot + "IOClockFlags",
                        network=network_["network"],
                        station=station_["station"],
                    )
                    tb2 = iof.to_dict("records")
                    cols2 = [{"name": i, "id": i} for i in iof.columns]

                else:
                    name = "Gap Metric Station Overview"
                    metric = db.get_metric(
                        metric_name=plot,
                        network=network_["network"],
                        station=station_["station"],
                    )
                    tb = metric.to_dict("records")
                    cols = [{"name": i, "id": i} for i in metric.columns]

            return (
                style_img,
                style_table,
                style_soh,
                style_soh,
                src,
                tb,
                cols,
                tb1,
                cols1,
                tb2,
                cols2,
                name,
                fig,
                style_graph,
            )


# ---------- Channel Functions -------------------

# This gets the channels available for specific network/station


@app.callback(
    Output("channel-selector", "options"),
    [
        Input("network", "children"),
        Input("station", "children"),
        Input("clicks", "children"),
    ],
    [State("input-db", "value")],
)
def _get_channel_selections(network, station, clicks, value):
    if clicks is not None and network is not None and station is not None:
        clicks_ = json.loads(clicks)
        network_ = json.loads(network)
        station_ = json.loads(station)
        if clicks_["n_clicks"] > 0:
            db = (
                Database(value)
                .view(network=network_["network"], station=station_["station"])
                .channel
            )
            db = db.where(db != "None")
            db = db.dropna()
            return [
                {"label": val, "value": val} for val in np.unique(db) if ":" not in val
            ]
        else:
            return []
    else:
        return []


# This stores the channel across tables


@app.callback(Output("channel", "children"), [Input("channel-selector", "value")])
def _get_channel(value):
    return json.dumps({"channel": value})


# ---- Channel Plot Selector Top Panel callbacks -----
## TODO: Need to add in location and time to display ##
# This gets the metrics available for a specific network/station/channel for the drop-downs in the top
# # panels of the channel page


@app.callback(
    Output("channel-plot-selectors-top-panel", "options"),
    [
        Input("network", "children"),
        Input("station", "children"),
        Input("channel", "children"),
        Input("clicks", "children"),
    ],
    [State("input-db", "value")],
)
def _get_metric_top_panel_selections(network, station, channel, clicks, value):
    if (
        clicks is not None
        and network is not None
        and station is not None
        and channel is not None
    ):
        clicks_ = json.loads(clicks)
        network_ = json.loads(network)
        station_ = json.loads(station)
        channel_ = json.loads(channel)
        if clicks_["n_clicks"] > 0:
            db = (
                Database(value)
                .view(
                    network=network_["network"],
                    station=station_["station"],
                    channel=channel_["channel"],
                )
                .metric
            )
            return [{"label": val, "value": val} for val in np.unique(db)]
        else:
            return [{"label": "", "value": ""}]
    else:
        return [{"label": "", "value": ""}]


# # Populate the table if the right metrics are selected, otherwise, populate a graph


@app.callback(
    [
        Output("top-metric-table", "data"),
        Output("top-metric-table", "columns"),
        Output("top-metric-table-div", "style"),
        Output("top-plot", "style_img"),
        Output("top-plot2", "style_img"),
        Output("div-image-holder-top", "style"),
        Output("top-plot", "src"),
        Output("top-plot2", "src"),
        Output("channel-top-graph", "figure"),
    ],
    [
        Input("network", "children"),
        Input("station", "children"),
        Input("channel", "children"),
        Input("clicks", "children"),
        Input("channel-plot-selectors-top-panel", "value"),
    ],
    [State("input-db", "value")],
)
# # Return appropriate metric table or plot based on metric selection
def _top_metric_table_plot(network, station, channel, clicks, metric, value):
    return _metric_table_plot(network, station, channel, clicks, metric, value)


# ---- Channel Plot Selector Botttom Panel callbacks -----
## TODO: Need to add in location and time to display ##

# This gets the metrics available for a specific network/station/channel for the drop-downs in the bottom
# panel of the channel page


@app.callback(
    Output("channel-plot-selectors-bottom-panel", "options"),
    [
        Input("network", "children"),
        Input("station", "children"),
        Input("channel", "children"),
        Input("clicks", "children"),
    ],
    [State("input-db", "value")],
)
def _get_metric_botttom_panel_selections(network, station, channel, clicks, value):
    if (
        clicks is not None
        and network is not None
        and station is not None
        and channel is not None
    ):
        clicks_ = json.loads(clicks)
        network_ = json.loads(network)
        station_ = json.loads(station)
        channel_ = json.loads(channel)
        if clicks_["n_clicks"] > 0:
            db = (
                Database(value)
                .view(
                    network=network_["network"],
                    station=station_["station"],
                    channel=channel_["channel"],
                )
                .metric
            )
            return [{"label": val, "value": val} for val in np.unique(db)]
        else:
            return [{"label": "", "value": ""}]
    else:
        return [{"label": "", "value": ""}]


# # Populate the table if the right metrics are selected, otherwise, populate a graph
@app.callback(
    [
        Output("bottom-metric-table", "data"),
        Output("bottom-metric-table", "columns"),
        Output("bottom-metric-table-div", "style"),
        Output("bottom-plot", "style"),
        Output("bottom-plot2", "style"),
        Output("div-image-holder-bottom", "style"),
        Output("bottom-plot", "src"),
        Output("bottom-plot2", "src"),
        Output("channel-bottom-graph", "figure"),
    ],
    [
        Input("network", "children"),
        Input("station", "children"),
        Input("channel", "children"),
        Input("clicks", "children"),
        Input("channel-plot-selectors-bottom-panel", "value"),
    ],
    [State("input-db", "value")],
)
# # Return appropriate metric table or plot based on metric selection
def _bottom_metric_table_plot(network, station, channel, clicks, metric, value):
    return _metric_table_plot(network, station, channel, clicks, metric, value)


def _metric_table_plot(network, station, channel, clicks, metric, value):
    if (
        clicks is not None
        and network is not None
        and station is not None
        and channel is not None
        and metric is not None
    ):
        clicks_ = json.loads(clicks)
        network_ = json.loads(network)
        station_ = json.loads(station)
        channel_ = json.loads(channel)
        metric_name = metric
        # Initialize values
        style = {"display": "none"}
        style_img = {"display": "none"}
        style_img2 = {"display": "none"}
        style_cont = {"display": "none"}
        columns = [
            {"name": i, "id": i}
            for i in [
                "network",
                "station",
                "channel",
                "location",
                "start_time",
                "end_time",
                "dc_offset_times",
                "masks",
            ]
        ]
        data = [
            {
                "network": "",
                "station": "",
                "channel": "",
                "location": "",
                "start_time": "",
                "end_time": "",
                "dc_offset_times": "",
                "masks": "",
            }
        ]
        src = ""
        src1 = ""
        fig = go.Figure()
        if clicks_["n_clicks"] > 0:
            if (
                metric_name == "deadChanMeanMetric"
                or metric_name == "deadChannelMetric"
            ):
                metric_name = "deadChannel"

            db = Database(value)
            net = network_["network"]
            sta = station_["station"]
            chan = channel_["channel"]
            if metric_name == "basicStatsMetric":
                style_img = {"display": "block"}
                style_cont = {"display": "block"}

                basic_stats_df = get_basic_stats_data(db, net, sta, chan)
                stat_list = [
                    "minimum",
                    "maximum",
                    "median",
                    "mean",
                    "variance",
                    "standard_deviation",
                    "rms",
                ]
                fig = make_subplots(
                    rows=7, cols=1, shared_xaxes=True, subplot_titles=tuple(stat_list)
                )
                row = 1
                basic_stats_df["start_time"] = basic_stats_df["start_time"].apply(
                    lambda x: x.date().strftime("%Y-%m-%d")
                )
                for stat in stat_list:
                    fig.append_trace(
                        go.Scatter(
                            x=basic_stats_df["start_time"],
                            y=basic_stats_df[stat],
                            name=stat,
                            mode="markers",
                        ),
                        row=row,
                        col=1,
                    )
                    row += 1
                fig.update_xaxes(
                    type="category",
                )
                fig.update_layout(height=800, width=600, showlegend=False)

            elif metric_name == "dailyPdfPlot":
                style_img = {"display": "block"}
                style_img2 = {"display": "block"}
                style_cont = {"display": "block"}

                df, df_line, snql = get_pdf_plot_data(
                    db, network=net, station=sta, channel=chan
                )
                df = plot_grid_data_fill_in(df)
                colormap = COLOR_SCHEMES["daily_pdf"]
                dates = df.axes[0].tolist()
                seconds = df.axes[1].tolist()
                x = [float(sec) for sec in seconds]
                y = dates
                z_vals = []
                for i in range(df.shape[0]):
                    col_vals = df.iloc[i]
                    z = []
                    for col_val in col_vals:
                        z.append(col_val)
                    z_vals.append(z)

                height_scale = 300 if len(y) == 1 else 600
                colorbar_scale = 2 if len(y) == 1 else 1
                fig = go.Figure(
                    data=go.Heatmap(
                        x=x,
                        y=y,
                        z=z_vals,
                        zmin=-10,
                        zmax=80,
                        zauto=False,
                        colorscale=colormap,
                        colorbar=dict(
                            len=colorbar_scale,
                            title="D",
                            titleside="top",
                            tickmode="array",
                            tickvals=np.linspace(-5, 75, 10),
                            ticktext=[
                                "D<0",
                                "10>D>=0",
                                "20>d>=10",
                                "30>D>=20",
                                "40>D>=30",
                                "50>D>=40",
                                "60>D>=50",
                                "70>D>=60",
                                "D>=70",
                                "No Data",
                            ],
                            ticks="outside",
                        ),
                    )
                )
                axis_template_y = dict(
                    showgrid=True,
                    zeroline=True,
                    linecolor="black",
                    showticklabels=True,
                    type="category",
                    ticks="outside",
                )

                axis_template_x = dict(
                    showgrid=True,
                    zeroline=True,
                    linecolor="black",
                    showticklabels=True,
                    ticks="outside",
                    tickvals=x,
                )

                fig.update_layout(
                    xaxis=axis_template_x,
                    yaxis=axis_template_y,
                    showlegend=True,
                    xaxis_title="Period (seconds)",
                    yaxis_title="Daily PDF Mode Power Grid: "
                    + net
                    + " "
                    + sta
                    + " "
                    + chan,
                    height=height_scale,
                    width=1600,
                )
                fig.update_xaxes(type="log")

            elif metric_name == "psdPlot":
                style_img = {"display": "block"}
                style_img2 = {"display": "block"}
                style_cont = {"display": "block"}
                (
                    psd,
                    psd_count,
                    period,
                    plot_style,
                    stats,
                    noise_matrix,
                    start_time,
                    end_time,
                    freq,
                ) = get_psd_plot_data(db, network=net, station=sta, channel=chan)
                psd_df = db.get_metric(
                    "psdPlot", network=net, station=sta, channel=chan
                )
                plot_style = str(psd_df["plot_style"])

                if "psd" in plot_style and len(psd) == 1:
                    fig = go.Figure()
                    for i in noise_matrix:
                        fig.add_trace(
                            go.Scatter(x=period, y=i, mode="lines", showlegend=False)
                        )
                    fig.update_xaxes(type="log")
                    fig.update_layout(
                        height=600,
                        width=1200,
                        title="{} Corrected, hourly PSDs from {} -- {} for {}.{}.{}".format(
                            psd_count, start_time, end_time, net, sta, chan
                        ),
                        xaxis_title="Period (s)",
                        yaxis_title="Power (dB)",
                    )

                elif "pdf" in plot_style:
                    pdf_matrix = np.fliplr(stats[0]["pdf_matrix"])
                    pdf_bins = stats[0]["pdf_bins"]
                    for indv, val in enumerate(pdf_matrix):
                        for inde, el in enumerate(val):
                            if el == 0.0:
                                pdf_matrix[indv][inde] = None
                    colormap = COLOR_SCHEMES["psd"]
                    fig = go.Figure(
                        data=go.Heatmap(
                            x=period[::-1],
                            y=pdf_bins,
                            z=pdf_matrix,
                            hoverongaps=False,
                            colorscale=colormap,
                        )
                    )
                    fig.update_xaxes(type="log", showgrid=True)
                    fig.update_layout(
                        height=600,
                        width=1200,
                        showlegend=True,
                        title="PDF plot of {} hourly PSDs from {} -- {} for {}.{}.{}".format(
                            psd_count, start_time, end_time, net, sta, chan
                        ),
                        xaxis_title="Period (s)",
                        yaxis_title="Amplitude [m^2/s^4/Hz] [dB]",
                    )

                if chan.startswith("B"):
                    fig.update_xaxes(range=(0.1, 100))
                freq_sub = freq[freq <= 10]
                period_sub = 1 / freq[freq <= 10]
                nlnm, nhnm = noiseModel(freq_sub)
                fig.add_trace(
                    go.Scatter(
                        x=period_sub,
                        y=nlnm,
                        name="NLNM",
                        line=dict(color="grey", width=2),
                    )
                )
                fig.add_trace(
                    go.Scatter(
                        x=period_sub,
                        y=nhnm,
                        name="NHNM",
                        line=dict(color="grey", width=2),
                    )
                )

                # This is logic for the default of showMedian being True
                fig.add_trace(
                    go.Scatter(
                        x=period,
                        # y=stats[0]["median"],
                        y=stats[0]["mode"],
                        # name="Median",
                        name="Mode",
                        line=dict(color="cyan", width=2),
                    )
                )

                fig.update_layout(
                    legend=dict(
                        orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1
                    )
                )

            else:
                style = {"display": "block"}

                # Connect to database and get
                db = Database(value).get_metric(
                    metric_name=metric_name,
                    network=network_["network"],
                    station=station_["station"],
                    channel=channel_["channel"],
                )
                # Create dataFrame
                dff = pd.DataFrame(db)
                # dff = db.drop(['created', 'metric_name'], axis=1)
                columns = [{"name": i, "id": i} for i in dff.columns]
                data = dff.to_dict("rows")

            return (
                data,
                columns,
                style,
                style_img,
                style_img2,
                style_cont,
                src,
                src1,
                fig,
            )
    else:
        style = {"display": "hidden"}
        style_img = {"display": "hidden"}
        style_img2 = {"display": "hidden"}
        style_cont = {"display": "hidden"}
        src = ""
        src1 = ""
        fig = go.Figure()
        columns = [
            {"name": i, "id": i}
            for i in [
                "network",
                "station",
                "channel",
                "location",
                "start_time",
                "end_time",
                "dc_offset_times",
                "masks",
            ]
        ]
        data = [
            {
                "network": "",
                "station": "",
                "channel": "",
                "location": "",
                "start_time": "",
                "end_time": "",
                "dc_offset_times": "",
                "masks": "",
            }
        ]
        return data, columns, style, style_img, style_img2, style_cont, src, src1, fig


if __name__ == "__main__":
    app.run_server(debug=True)
