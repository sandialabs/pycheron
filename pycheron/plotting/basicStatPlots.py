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

import pandas as pd
import os
import matplotlib.pyplot as plt
import csv
from pycheron.db.sqllite_db import Database

__all__ = ["basicstat_timeseries_plot", "basicstat_box_plot", "get_basic_stats_data"]

y_vals = [
    "rms",
    "standard_deviation",
    "median",
    "maximum",
    "minimum",
    "variance",
    "mean",
]


def basicstat_timeseries_plot(data, output_dir):
    """
    Creates a time series plot of basic stats (1 plot per metric) for each station-channel pair

    :param data: Database object containing metrics to plot
    :type data: Pycheron.db.sqllite_db.Database
    :param output_dir: Directory in which the plots will be saved
    :type output_dir: str

    :return: returns time series plot of basic stats (1 plot per metric/station chan) over a given time period within
             the input database

    """
    if isinstance(data, Database):
        networks = data.networks()
        stations = data.stations()
        channels = data.channels()
        for network in networks:
            for station in stations:
                for channel in channels:
                    basic_stat_df = get_basic_stats_data(
                        data, network, station, channel
                    )
                    _make_time_series_plot(
                        basic_stat_df, output_dir, network, station, channel
                    )

    return


def _make_time_series_plot(basic_stat_df, output_dir, network, station, channel):
    title = "{}.{}.{} basic stats".format(network, station, channel)
    basic_stat_df.plot(
        x="start_time",
        y=y_vals,
        style=".",
        subplots=True,
        grid=True,
        title=title,
        layout=(7, 1),
        sharex=True,
        sharey=False,
        legend=True,
    )

    # Save figure to a png
    png_name = (
        output_dir + "/" + "{}.{}.{}scatter.png".format(network, station, channel)
    )
    plt.savefig(png_name, dpi=500)
    plt.close()


def get_basic_stats_data(data, network, station, channel):
    basic_stat_df = data.get_metric(
        "basicStatsMetric",
        network=network,
        station=station,
        channel=channel,
    )
    if not basic_stat_df.empty:
        for y in y_vals:
            basic_stat_df[y] = pd.to_numeric(basic_stat_df[y])
            # Convert start_time column to datetime pandas can plot
        basic_stat_df["start_time"] = pd.to_datetime(basic_stat_df["start_time"])

    return basic_stat_df


def basicstat_box_plot(data, output_dir):
    """
    Creates a box plot of basic stats (1 plot per metric) for each channel

    :param data: Database object containing metrics to plot
    :type data: Pycheron.db.sqllite_db.Database
    :param output_dir: Directory in which the plots will be saved
    :type output_dir: str
    :return: returns box plot of basic stats (1 plot per metric/chan) over a given time period within the input
             database
    """

    if isinstance(data, Database):
        networks = data.networks()
        stations = data.stations()
        channels = data.channels()
        for network in networks:
            for station in stations:
                for channel in channels:
                    basic_stat_df = get_basic_stats_data(
                        data, network, station, channel
                    )
                    _make_box_plots(
                        basic_stat_df, output_dir, network, station, channel
                    )

    return


def _make_box_plots(basic_stat_df, output_dir, network, station, channel):
    # Plot subplot of 1,7 with each boxplot in one pane
    fig, axes = plt.subplots(7, 1, figsize=(10, 12))
    fig.tight_layout()
    basic_stat_df.boxplot(column="rms", ax=axes[0])
    basic_stat_df.boxplot(column="standard_deviation", ax=axes[1])
    basic_stat_df.boxplot(column="median", ax=axes[2])
    basic_stat_df.boxplot(column="maximum", ax=axes[3])
    basic_stat_df.boxplot(column="minimum", ax=axes[4])
    basic_stat_df.boxplot(column="variance", ax=axes[5])
    basic_stat_df.boxplot(column="mean", ax=axes[6])

    # Save figure to a png
    png_name = (
        output_dir + "/" + "{}.{}.{}boxplot.png".format(network, station, channel)
    )
    plt.savefig(png_name, dpi=500)
    plt.close()
