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

import json

import numpy as np
import pandas as pd

from pycheron.db.sqllite_db import Database
from pycheron.plotting.dailyPdfPlot import noiseDiff, plot_grid, _find_nearest
from pycheron.psd.noise.networkNoiseModel import networkNoiseModel
from pycheron.psd.psdList import psdList
from pycheron.psd.psdStatistics import psdStatistics
from pycheron.util.format import parse_snclq
from pycheron.util.getChannel import getChannelName
from pycheron.util.getMPeriod import get_m_period
from pycheron.util.getSta import getStationName

from pycheron.plotting.psdPlot import _read_psds

__all__ = ["rankPlot", "get_rank_day_averages"]


def rankPlot(
    st,
    f_name=None,
    model="nlnm",
    rank_by=None,
    nnm_fname=None,
    database=None,
    station=None,
    network=None,
    channel=None,
    location=None,
    session=None,
):
    """
    This takes in network level streams and ranks the plots

    :param st: ObsPy stream OR Database object
    :type st: obspy.core.stream.Stream or pycheron.db.sqllite_db.Database
    :param f_name: name of output file to save ranking plot
    :type f_name: str
    :param model: Model. "nlnm" (new low noise model) or "nnm" (network noise model).
    :type model: str
    :param rank_by: Array of period values to rank the stations at. If none, default will be [1, 6.5, 30, 100]
    :type rank_by: numpy.ndarray
    :param nnm_fname: If supplied, plot of network noise model will be saved
    :type nnm_fname: str
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database
    :param station: Use only if **st** is database object. Station name
    :type station: str
    :param network: Use only if **st** is database object. Network name
    :type network: str
    :param channel: Use only if **st** is database object. Channel name
    :type channel: str
    :param location: Use only if **st** is database object. Location name
    :type location: str
    :param session: Use only if **st** is database object. Session name
    :type session: str

    Station ranking plots were first established for use in the
    USArray TA project (Busby et. al., 2018) and code was ported from Perl and Hypertext Preprocessor (PHP) to Python.
    Code was adapted and augmented for use within Pycheron.

    #. Busby, R. W., R.L. Woodward, K.A. Hafner, F. L. Vernon, and A.M. Frasseto (2018). The Design and Implementation
    of EarthScope's USArray Transportable Array in the Conterminous United States and Southern Canada (Rep.).


    **Example**

    .. code-block:: python

        # Import libs
        from obspy.clients.fdsn import Client
        from plotting.rankingPlot import rankPlot
        from obspy import UTCDateTime

        client = Client("IRIS")
        t = UTCDateTime("2019-09-10T00:00:00")
        te = UTCDateTime("2019-09-10T05:00:00")
        st = client.get_waveforms("CI", "*", "*", "H*", t, te)

        rankPlot(st, model="gsn")

    .. image:: _static/rankPlot.png

    """

    # If rank not provided by users
    if rank_by is None:
        rank_by = [1, 6.5, 30, 100]

    # Create dataframe for each component
    df_e = pd.DataFrame()
    df_n = pd.DataFrame()
    df_z = pd.DataFrame()

    # If database instance call with appropriate args
    if isinstance(st, Database):
        _rankPlot_from_database(
            st,
            f_name,
            model,
            rank_by,
            nnm_fname,
            station,
            network,
            channel,
            location,
            session,
        )
    # Otherwise parse snclq information
    else:
        network, station, channel, location, quality = parse_snclq(st[0].get_id())
        snclq = st[0].get_id()
        # Network Noise Model
        if model == "nnm":
            # If database object, check to see if we already have the psds calculated
            if database is not None:
                time = (
                    st[0].stats.starttime.isoformat(),
                    st[0].stats.endtime.isoformat(),
                )
                tb = database.get_metric("networkNoiseModel", network=network, session=session)

                # if returns no data, calculate and insert into db
                if tb.empty:
                    networkNoiseModel(
                        st,
                        plot=True,
                        fname=nnm_fname,
                        station=station,
                        network=network,
                        channel=channel,
                        session=session,
                        database=database,
                    )

                    _rankPlot_from_database(
                        database,
                        f_name,
                        model,
                        rank_by,
                        nnm_fname,
                        station,
                        network,
                        channel,
                        location,
                        session,
                    )
                else:
                    _rankPlot_from_database(
                        database,
                        f_name,
                        model,
                        rank_by,
                        nnm_fname,
                        station,
                        network,
                        channel,
                        location,
                        session,
                    )
            # Otherwise, calculate a network noise model
            else:
                nnm_model = networkNoiseModel(st, plot=True, fname=nnm_fname)
                # Stations for E/N Model
                sta_en = nnm_model["stations_en"]
                # Stations for Z Model
                sta_z = nnm_model["stations_z"]
                # Power from Z, E/N Model
                m_power_z = nnm_model["zModel"][1]
                m_power_en = nnm_model["enModel"][1]
                # Get period values
                m_period = get_m_period(nnm_model)

                # Grabbing psd stats from output of networkNoiseModel
                psds = nnm_model["psdStats"]  # TODO
                # PSD periods
                period = 1 / psds["noise_matrix_frequency"][0][0][np.where(m_period)]

                # for Z channel
                # looping through channels
                for i in range(len(psds["snclq"][0])):
                    # loopinf through Z stations
                    for j in range(len(sta_z)):
                        # grabbing channel name
                        chan = getChannelName(psds["snclq"][j][i])
                        # grabbing Z channel
                        if chan.endswith("Z"):
                            chan_z = chan
                            # getting mode
                            mode = psds["mode"][j][i]
                            # getting the difference between actual and network noise model
                            diff, per = noiseDiff(period, mode, m_power_z, m_period)
                            # adding data to DF
                            df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                            df_z = pd.concat([df_z, df2], ignore_index=True)

                # adding indexs
                df_z["Station"] = sta_z
                df_z = df_z.rename(df_z["Station"])
                df_z = df_z.drop(["Station"], axis=1)

                for i in range(len(psds["snclq"][0])):
                    # looping through E/N stations
                    for j in range(len(sta_en)):
                        # grabbing channel name
                        chan = getChannelName(psds["snclq"][j][i])
                        # grabbing E channel
                        if chan.endswith("E"):
                            chan_e = chan
                            mode = psds["mode"][j][i]
                            # getting the difference between actual and network noise model
                            diff, per = noiseDiff(period, mode, m_power_en, m_period)
                            # adding data to DF
                            df2_e = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                            df_e = pd.concat([df_e, df2_e], ignore_index=True)

                        # grabbing N channel
                        if chan.endswith("N"):
                            chan_n = chan
                            mode = psds["mode"][j][i]
                            # getting the difference between actual and network noise model
                            diff, per = noiseDiff(period, mode, m_power_en, m_period)
                            # adding data to DF
                            df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                            df_n = pd.concat([df_n, df2], ignore_index=True)
                # adding station as index
                df_e["Station"] = sta_en
                df_e = df_e.rename(df_e["Station"])
                df_e = df_e.drop(["Station"], axis=1)
                # adding station as index
                df_n["Station"] = sta_en
                df_n = df_n.rename(df_n["Station"])
                df_n = df_n.drop(["Station"], axis=1)
        # If NLNM model
        if model == "nlnm":
            # Get time
            if database is not None:
                time = (
                    st[0].stats.starttime.isoformat(),
                    st[0].stats.endtime.isoformat(),
                )
                # Get metric information if database available
                tb = database.get_metric(
                    "networkNoiseModel",
                    network=network,
                    station=station,
                    channel=channel,
                    location=location,
                    session=session,
                    time=time,
                )

                # if returns no data, calculate and insert into db
                if tb.empty:
                    networkNoiseModel(
                        st,
                        plot=True,
                        fname=nnm_fname,
                        station=station,
                        network=network,
                        channel=channel,
                        session=session,
                        database=database,
                    )

                    _rankPlot_from_database(
                        st,
                        f_name,
                        model,
                        rank_by,
                        nnm_fname,
                        station,
                        network,
                        channel,
                        location,
                        session,
                    )
                # Otherwise just get data
                else:
                    _rankPlot_from_database(
                        st,
                        f_name,
                        model,
                        rank_by,
                        nnm_fname,
                        station,
                        network,
                        channel,
                        location,
                        session,
                    )
            # If not database, then start calculating everything
            else:
                # Initialize lists for each component
                sta_n = []
                sta_e = []
                sta_z = []
                # Calculate psds, get freq, period
                psds = psdList(st)
                freq = psds[0][0][0]
                period = 1 / freq
                period = period[np.where(period > 0.1)]
                # Calculate psd statistics
                stats = psdStatistics(psds, evalresp=None, database=database)
                # Loop through the stats and grab out the sta, chan, mode, period and power and difference between NLNM
                # for each component
                for i in range(len(stats)):
                    chan = getChannelName(stats[i]["snclq"])
                    # N component
                    if chan.endswith("N"):
                        sta_n.append(stats[i]["snclq"].split(".")[1])
                        chan_n = chan
                        mode = stats[i]["mode"][np.where(period)]
                        m_period = 1 / stats[i]["noise_matrix_frequency"][np.where(period)]
                        m_power = stats[i]["nlnm"][np.where(period)]
                        # Get diff between model and psds
                        diff, per = noiseDiff(period, mode, m_power, m_period)
                        # adding data to DF
                        df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                        df_n = pd.concat([df_n, df2], ignore_index=True)
                    # E component
                    if chan.endswith("E"):
                        sta_e.append(stats[i]["snclq"].split(".")[1])
                        chan_e = chan
                        mode = stats[i]["mode"][np.where(period)]
                        m_period = 1 / stats[i]["noise_matrix_frequency"][np.where(period)]
                        m_power = stats[i]["nlnm"][np.where(period)]
                        diff, per = noiseDiff(period, mode, m_power, m_period)
                        # adding data to DF
                        df2_e = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                        df_e = pd.concat([df_e, df2_e], ignore_index=True)
                    # Z component
                    if chan.endswith("Z"):
                        sta_z.append(stats[i]["snclq"].split(".")[1])
                        chan_z = chan
                        mode = stats[i]["mode"][np.where(period)]
                        m_period = 1 / stats[i]["noise_matrix_frequency"][np.where(period)]
                        m_power = stats[i]["nlnm"][np.where(period)]
                        diff, per = noiseDiff(period, mode, m_power, m_period)
                        # adding data to DF
                        df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                        df_z = pd.concat([df_z, df2], ignore_index=True)

                # adding station as index
                df_z["Station"] = sta_z
                df_z = df_z.rename(df_z["Station"])
                df_z = df_z.drop(["Station"], axis=1)
                # adding station as index
                df_e["Station"] = sta_e
                df_e = df_e.rename(df_e["Station"])
                df_e = df_e.drop(["Station"], axis=1)
                # adding station as index
                df_n["Station"] = sta_n
                df_n = df_n.rename(df_n["Station"])
                df_n = df_n.drop(["Station"], axis=1)

        # ranking data
        for k in rank_by:
            # finding period closest to input
            if not df_z.empty:
                # finding period closest to input
                index_dz = _find_nearest(np.asarray(df_z.columns, dtype=float), k)
                # sorting values
                df_z_rank = df_z.sort_values(by=[df_z.columns[index_dz]])
                # plotting
                if f_name is None:
                    f_name = "f_name"
                # Create plot
                plot_grid(
                    df_z_rank,
                    snclq,
                    plot="stationRanking",
                    f_name=f_name + chan_z + "_" + str(k),
                    channel=chan_z,
                    rank=float(df_z.columns[index_dz]),
                )

            if not df_e.empty:
                # finding period closest to input
                index_de = _find_nearest(np.asarray(df_e.columns, dtype=float), k)
                # sorting values
                df_e_rank = df_e.sort_values(by=[df_e.columns[index_de]])
                # plotting
                if f_name is None:
                    f_name = "f_name"
                # Create plot
                plot_grid(
                    df_e_rank,
                    snclq,
                    plot="stationRanking",
                    f_name=f_name + chan_e + "_" + str(k),
                    channel=chan_e,
                    rank=float(df_e.columns[index_de]),
                )

            if not df_n.empty:
                # finding period closest to input
                index_dn = _find_nearest(np.asarray(df_n.columns, dtype=float), k)
                # sorting values
                df_n_rank = df_n.sort_values(by=[df_n.columns[index_dn]])
                # plotting
                if f_name is None:
                    f_name = "f_name"
                # Create plot
                plot_grid(
                    df_n_rank,
                    snclq,
                    plot="stationRanking",
                    f_name=f_name + chan_n + "_" + str(k),
                    channel=chan_n,
                    rank=float(df_n.columns[index_dn]),
                )


def _rankPlot_from_database(
    database,
    f_name=None,
    model="nlnm",
    rank_by=None,
    nnm_fname=None,
    station=None,
    network=None,
    channel=None,
    location=None,
    session=None,
):
    """
    Internal function to generate ranking plot from existing database information
    """
    df_e = pd.DataFrame()
    df_n = pd.DataFrame()
    df_z = pd.DataFrame()

    if model == "nnm":

        tb = database.get_metric("networkNoiseModel", network=network, session=session)

        # if returns no data, calculate and insert into db
        if tb.empty:
            networkNoiseModel(
                database,
                plot=True,
                fname=nnm_fname,
                station=station,
                network=network,
                channel=channel,
                session=session,
                database=database,
            )

            tb = database.get_metric("networkNoiseModel", network=network, session=session)

        # find index number of nnm which has highest number of stations
        num_stations_en = 0
        index_en = 0
        num_stations_z = 0
        index_z = 0
        for i in range(len(tb)):
            stations_en = len(np.unique(json.loads(tb["stations_en"].iloc[i])))
            stations_z = len(np.unique(json.loads(tb["stations_z"].iloc[i])))
            if num_stations_en < stations_en:
                num_stations_en = stations_en
                index_en = i
            if num_stations_z < stations_z:
                num_stations_z = stations_z
                index_z = i
        # Stations for E/N Model
        sta_en = json.loads(tb["stations_en"].iloc[index_en])
        # Stations for Z Model
        sta_z = json.loads(tb["stations_z"].iloc[index_z])
        # Power from E/N Model

        m_power_en = json.loads(tb["enModel"].iloc[index_en])[1]
        # Power from Z Model
        m_power_z = json.loads(tb["zModel"].iloc[index_z])[1]
        # Model period
        # add if to determine between z period and en model period
        if num_stations_z > 0 and num_stations_en == 0:
            m_period = json.loads(tb["zModel"].iloc[index_z])[0]
            stats = json.loads(tb["psdStats"].iloc[index_z])
        else:
            m_period = json.loads(tb["enModel"].iloc[index_en])[0]
            stats = json.loads(tb["psdStats"].iloc[index_en])
        for i in range(len(stats["mode"])):
            stats["snclq"][i] = json.loads(stats["snclq"][i])
            for j in range(len(stats["mode"][i])):
                stats["noise_matrix_frequency"][i][j] = json.loads(stats["noise_matrix_frequency"][i][j])
                stats["mode"][i][j] = json.loads(stats["mode"][i][j])

        # PSD periods
        freq = np.asarray(stats["noise_matrix_frequency"][0][0])
        mode_list = stats["mode"]
        snclq_list = stats["snclq"]
        snclq = stats["snclq"][0]  # need for color grad but only uses network, so does not matter which snclq we use
        period = 1 / freq[np.where(m_period)]
        # for Z channel
        # looping through stations in snclq list
        for i in range(len(sta_z)):
            # looping through each channel
            df2 = pd.DataFrame()
            for j in range(len(snclq_list)):
                if sta_z[i] == getStationName(snclq_list[j][0]):
                    # grabbing channel name
                    for k in range(len(snclq_list[j])):
                        chan = getChannelName(snclq_list[j][k])
                        # grabbing Z channel
                        if chan.endswith("Z"):
                            chan_z = chan
                            # getting mode
                            mode = mode_list[j][k]
                            # getting the difference between actual and network noise model
                            diff, per = noiseDiff(period, mode, m_power_z, m_period)
                            # adding data to DF
                            df_z_temp = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                            df2 = pd.concat([df2, df_z_temp])
            df2 = df2.mode()
            if len(df2) > 1:
                df2 = df2.mean()
            df_z = pd.concat([df_z, df2], ignore_index=True)
        # adding indexs
        df_z["Station"] = sta_z
        df_z = df_z.rename(df_z["Station"])
        df_z = df_z.drop(["Station"], axis=1)

        for i in range(len(sta_en)):
            # loopinf through E/N stations
            df2_e = pd.DataFrame()
            df2_n = pd.DataFrame()
            for j in range(len(snclq_list)):
                if sta_en[i] == getStationName(snclq_list[j][0]):
                    for k in range(len(snclq_list[j])):
                        # grabbing channel name
                        chan = getChannelName(snclq_list[j][k])
                        # grabbing E channel
                        if chan.endswith("E"):
                            chan_e = chan
                            mode = mode_list[j][k]
                            # getting the difference between actual and network noise model
                            diff, per = noiseDiff(period, mode, m_power_en, m_period)
                            # adding data to DF
                            df2_e_temp = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                            df2_e = pd.concat([df2_e, df2_e_temp])

                        # grabbing N channel
                        if chan.endswith("N"):
                            chan_n = chan
                            mode = mode_list[j][k]
                            # getting the difference between actual and network noise model
                            diff, per = noiseDiff(period, mode, m_power_en, m_period)
                            # adding data to DF
                            df2_n_temp = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                            df2_n = pd.concat([df2_e, df2_n_temp])

            df2_e = df2_e.mode()
            if len(df2_e) > 1:
                df2_e = df2_e.mean()
            df_e = pd.concat([df_e, df2_e], ignore_index=True)

            df2_n = df2_n.mode()
            if len(df2_n) > 1:
                df2_n = df2_n.mean()
            df_n = pd.concat([df_n, df2_n], ignore_index=True)

        # adding station as index
        df_e["Station"] = sta_en
        df_e = df_e.rename(df_e["Station"])
        df_e = df_e.drop(["Station"], axis=1)
        # adding station as index
        df_n["Station"] = sta_en
        df_n = df_n.rename(df_n["Station"])
        df_n = df_n.drop(["Station"], axis=1)

    if model == "nlnm":
        psd_channels, stats, period = calc_stats_from_psds_rankplot(
            database, network, channel, station, location, session
        )
        df_dict = calc_power_period_rankplot(stats, period)
        df_z, df_e, df_n = drop_and_rename_rankplot(df_dict)
        chan_z, chan_e, chan_n = df_dict["chan_z"], df_dict["chan_e"], df_dict["chan_n"]
        snclq = str(psd_channels[0][0][0][0][0][2][0])

    # ranking data
    for k in rank_by:
        # finding period closest to input
        if df_z.empty is False:
            # finding period closest to input
            index_dz = _find_nearest(np.asarray(df_z.columns, dtype=float), k)
            # sorting values
            df_z_rank = df_z.sort_values(by=[df_z.columns[index_dz]])
            # plotting
            if f_name is None:
                f_name = "f_name"
            plot_grid(
                df_z_rank,
                snclq,
                plot="stationRanking",
                f_name=f_name + chan_z + "_" + str(k),
                channel=chan_z,
                rank=float(df_z.columns[index_dz]),
            )

        if df_e.empty is False:
            # finding period closest to input
            index_de = _find_nearest(np.asarray(df_e.columns, dtype=float), k)
            # sorting values
            df_e_rank = df_e.sort_values(by=[df_e.columns[index_de]])
            # plotting
            if f_name is None:
                f_name = "f_name"
            plot_grid(
                df_e_rank,
                snclq,
                plot="stationRanking",
                f_name=f_name + chan_e + "_" + str(k),
                channel=chan_e,
                rank=float(df_e.columns[index_de]),
            )

        if df_n.empty is False:
            # finding period closest to input
            index_dn = _find_nearest(np.asarray(df_n.columns, dtype=float), k)
            # sorting values
            df_n_rank = df_n.sort_values(by=[df_n.columns[index_dn]])
            # plotting
            if f_name is None:
                f_name = "f_name"
            plot_grid(
                df_n_rank,
                snclq,
                plot="stationRanking",
                f_name=f_name + chan_n + "_" + str(k),
                channel=chan_n,
                rank=float(df_n.columns[index_dn]),
            )


def calc_stats_from_psds_rankplot(database, network, channel, station=None, location=None, session=None):
    """
    Function to calculate psd statistics from psds that exist within the database
    """

    # Grab out psds from database
    psd_channels = _read_psds(
        database,
        network=network,
        station=station,
        channel=channel,
        location=location,
        session=session,
    )

    # If psd_channels is empty exit
    if not psd_channels:
        return

    # Get freq, period
    freq = np.array(psd_channels[0][0][0][0][0][0])
    period = 1 / freq
    period = period[np.where(period > 0.1)]
    stats = psdStatistics(psd_channels[0][0][0], evalresp=None)
    return (psd_channels, stats, period)


def calc_power_period_rankplot(stats, period):
    """
    Function to calculate power and period for rank plots given psd statistics and period
    """

    # Create empty station, channel component lists and dataframes
    sta_n, sta_e, sta_z = [], [], []
    chan_n, chan_e, chan_z = None, None, None
    df_n, df_e, df_z = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    # Loop through psd statistics for each channel
    for i in range(len(stats)):
        chan = getChannelName(stats[i]["snclq"])
        stats[i]["noise_matrix_frequency"] = np.array(stats[i]["noise_matrix_frequency"])
        # N component
        if chan.endswith("N"):
            # Get sta, channel
            sta_n.append(stats[i]["snclq"].split(".")[1])
            chan_n = chan
            # Get mode, model period and model power
            mode = stats[i]["mode"][np.where(period)]
            m_period = 1 / stats[i]["noise_matrix_frequency"][np.where(period)]
            m_power = stats[i]["nlnm"][np.where(period)]
            # Difference psds from model
            diff, per = noiseDiff(period, mode, m_power, m_period)
            # adding data to DF
            df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
            df_n = pd.concat([df_n, df2], ignore_index=True)
        # E component
        if chan.endswith("E"):
            # Get sta, channel
            sta_e.append(stats[i]["snclq"].split(".")[1])
            chan_e = chan
            # Get mode, model period and model power
            mode = stats[i]["mode"][np.where(period)]
            m_period = 1 / stats[i]["noise_matrix_frequency"][np.where(period)]
            m_power = stats[i]["nlnm"][np.where(period)]
            # Difference psds from model
            diff, per = noiseDiff(period, mode, m_power, m_period)
            # adding data to DF
            df2_e = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
            df_e = pd.concat([df_e, df2_e], ignore_index=True)
        # Z component
        if chan.endswith("Z"):
            # Get sta, channel
            sta_z.append(stats[i]["snclq"].split(".")[1])
            chan_z = chan
            # Get mode, model period and model power
            mode = stats[i]["mode"][np.where(period)]
            m_period = 1 / stats[i]["noise_matrix_frequency"][np.where(period)]
            m_power = stats[i]["nlnm"][np.where(period)]
            # Difference psds from model
            diff, per = noiseDiff(period, mode, m_power, m_period)
            # adding data to DF
            df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
            df_z = pd.concat([df_z, df2], ignore_index=True)
    df_dict = {}
    df_dict["sta_n"], df_dict["sta_e"], df_dict["sta_z"] = sta_n, sta_e, sta_z
    df_dict["chan_n"], df_dict["chan_e"], df_dict["chan_z"] = chan_n, chan_e, chan_z
    df_dict["df_n"], df_dict["df_e"], df_dict["df_z"] = df_n, df_e, df_z

    return df_dict


def drop_and_rename_rankplot(df_dict):
    """
    Simple function to add station as an index to each channel component dataframe
    """
    # adding station as index
    df_dict["df_z"]["Station"] = df_dict["sta_z"]
    df_dict["df_z"] = df_dict["df_z"].rename(df_dict["df_z"]["Station"])
    df_dict["df_z"] = df_dict["df_z"].drop(["Station"], axis=1)
    # adding station as index
    df_dict["df_e"]["Station"] = df_dict["sta_e"]
    df_dict["df_e"] = df_dict["df_e"].rename(df_dict["df_e"]["Station"])
    df_dict["df_e"] = df_dict["df_e"].drop(["Station"], axis=1)
    # adding station as index
    df_dict["df_n"]["Station"] = df_dict["sta_n"]
    df_dict["df_n"] = df_dict["df_n"].rename(df_dict["df_n"]["Station"])
    df_dict["df_n"] = df_dict["df_n"].drop(["Station"], axis=1)
    return df_dict["df_z"], df_dict["df_e"], df_dict["df_n"]


def get_rank_day_averages(df_rank):
    """
    Simple function to get the mean value of each day for ranking plots
    """
    # Get mean values of days
    unique_frame = pd.DataFrame(columns=df_rank.columns, index=df_rank.index.unique())
    for name in df_rank.index.unique():
        if df_rank.loc[name].mean().size > 1:
            unique_frame.loc[name] = df_rank.loc[name].mean().values
        else:
            unique_frame.loc[name] = df_rank.loc[name]
    return unique_frame
