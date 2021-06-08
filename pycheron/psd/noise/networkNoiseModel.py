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

import matplotlib.pyplot as plt
import numpy as np
import obspy
import pandas as pd
from collections import OrderedDict

from obspy import UTCDateTime
from obspy.core.stream import Stream

from pycheron.db.sqllite_db import Database
from pycheron.psd.psdList import psdList
from pycheron.psd.psdStatistics import psdStatistics

from pycheron.util.logger import Logger
from pycheron.util.getSta import getSta, getUniqSta
from pycheron.util.getChannel import getChannelName
from pycheron.util.format import parse_snclq

plt.style.use("ggplot")
__all__ = ["networkNoiseModel"]


def networkNoiseModel(
    st,
    plot=False,
    fname=None,
    evalresp=None,
    station=None,
    network=None,
    channel=None,
    location=None,
    time=None,
    ylo=None,
    yhi=None,
    session=None,
    database=None,
    logger=None,
):
    """
    Calculates Network Noise Model for network.

    :param st: Stream of data that includes stations from network. Recommended that length > 1 Week
                OR Database object
    :type st: obspy.core.stream.Stream or pycheron.db.sqllite_db.Database
    :param plot: To plot or not (Default = False)
    :type plot: bool
    :param fname: If plotting at want to save the output
    :type fname: bool
    :param evalresp: IRIS evalresp file or array
    :type evalresp: numpy.ndarray
    :param network: If using database, **must** supply the network
    :type network: str
    :param station: If using database, and want to use specific stations for noise model
    :type station: str
    :param channel: If using database, and want to use specific channels for noise model
    :type channel: str
    :param session: specific session in database
    :type session: str
    :param database: Database object to save to
    :type database: pycheron.db.sqllite_db.Database
    :param logger: Logger object
    :type logger: pycheron.util.logger.Logger
    :param location: If using database, and want to use specific location for noise model
    :type location: str

    :return:

        * Dataframe of Network Noise Model from all \*\*E/\*\*N Stations
        * Dataframe of Network Noise Model from all \*\*Z Stations

    :rtype: dict

    Network Noise Model concept was borrowed from the USArray TA project (Busby et. al., 2018)

    #. Busby, R. W., R.L. Woodward, K.A. Hafner, F. L. Vernon, and A.M. Frasseto (2018). The Design and Implementation
    of EarthScope's USArray Transportable Array in the Conterminous United States and Southern Canada (Rep.).

    **Algorithm Steps of Network Noise Model**

    1) Calculate the mode at all periods for each station for extent of data provided
    2) Remove stations that have any mode values less than the NLNM (e.g., returning digitizer noise),
       mode values > -80db
    3) Take the mean of the modes at each of the periods.
    4) E/N channels are combined into a single horizontal model; Z is its own vertical model

    **Example**

    .. code-block:: python

        from pycheron.psd.noise.networkNoiseModel import networkNoiseModel
        from obspy import UTCDateTime
        from obspy.clients.fdsn import Client

        client = Client("IRIS")
        t = UTCDateTime("2017-09-10T00:00:00")
        te = UTCDateTime("2017-09-11T00:00:00")
        st = client.get_waveforms("AG", "*", "*", "*", t, te)

        model = networkNoiseModel(st)

    **Plotting**

    .. code-block:: python

        import matplotlib.pyplot as plt

        # uses output from above example
        plt.semilogx(model["enModel"], 'orangered', label="E/N Network Model")

        plt.semilogx(model["zModel"], 'c', label="Z Network Model")

        plt.xlabel("Period")
        plt.ylabel("Power(dB)")
        plt.title("Noise Model for Network: " + st[0].get_id().split(".")[0])
        plt.grid(True, "both", "both")
        plt.legend()

    .. image:: _static/networkNoiseModel.png
    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Initiate empty pandas dataframes
    df_sta = pd.DataFrame()
    dfz_sta = pd.DataFrame()

    # If dealing with obspy stream, get unique stations and network name, then loop through stations
    if isinstance(st, obspy.core.stream.Stream):
        st = multiday_traces_to_single_days(st)
        sta = getUniqSta(st)
        network = st[0].get_id().split(".")[0]
        # Loop through unique stations
        for i in range(len(sta)):
            # Obtain all traces in stream that match the specified station
            stN = getSta(st, sta[i])
            # parse snclq
            network, station, channel, location, quality = parse_snclq(stN[0].get_id())
            # If database object, check to see if we already have the psds calculated
            if database is not None:
                time = (
                    st[0].stats.starttime.isoformat(),
                    st[0].stats.endtime.isoformat(),
                )
                tb = database.get_metric(
                    "psdMetric",
                    network=network,
                    station=station,
                    channel=channel,
                    location=location,
                    session=session,
                    time=time,
                )

                # if returns no data, calculate psds and insert into db
                if tb.empty:
                    psds = psdList(stN)
                    database.insert_metric(psds)

                else:
                    tb_psds = tb.uncorrected_psds[tb.channel == channel]
                    psds = []
                    for j in range(len(tb_psds)):
                        psd = database.extract_masks(tb_psds.iloc[j])
                        psds.append(psd)
            # If not database object, calculate psds and psd Statistics
            else:
                psds = psdList(stN)

            stats = psdStatistics(psds, logger=logger, evalresp=evalresp, database=database)
            # get period where freq < 10Hz
            try:
                freq = np.asarray(psds[0][0][0])
                per = 1 / freq
            except IndexError:
                # del sta_n[i]
                continue
            # get NLNM where freq < 10Hz
            nlnm = stats[0]["nlnm"][np.where(per)]
            # get channels for df  index
            chan = []
            # extract each channel modes
            df_chan = pd.DataFrame()
            dfz_chan = pd.DataFrame()
            for j in range(len(stats)):
                # getting mode for each channel
                try:
                    mode = stats[j]["mode"][np.where(per)]
                except IndexError:
                    continue
                chan.append(getChannelName(stats[j]["snclq"]))
                # adding modes to df
                df2 = pd.DataFrame(data=mode.reshape(1, len(mode)), columns=per)
                df_chan = df_chan.append(df2, ignore_index=True)
            # adding indexes
            df_chan["snclq"] = chan

            # creating copy of df, so wont get index error as I am dropping Z channels
            df_chan_copy = df_chan
            # loop through to find z channel

            for k in range(len(df_chan.index)):
                name = df_chan.iloc[k].snclq
                if name.endswith("Z"):
                    # make z channel have own df
                    df = pd.DataFrame(
                        df_chan.iloc[k].values.reshape(1, len(df_chan.iloc[k].values)),
                        columns=list(df_chan.iloc[k].index),
                    )
                    dfz_chan = dfz_chan.append(df, ignore_index=True)
                    # drop z channel from E/N df
                    df_chan_copy = df_chan_copy.drop(k)
            # reassigning the df with dropped Z channels to original df
            df_chan = df_chan_copy

            df_chan = df_chan.rename(df_chan["snclq"])
            df_chan = df_chan.drop(["snclq"], axis=1)

            if not df_chan.empty:
                # getting the mean between E/N station modes
                mode_sta_en = df_chan.mean()
            else:
                mode_sta_en = pd.DataFrame()

            if not dfz_chan.empty:
                dfz_chan = dfz_chan.rename(dfz_chan["snclq"])
                dfz_chan = dfz_chan.drop(["snclq"], axis=1)
                mode_sta_z = dfz_chan.mean()
            else:
                mode_sta_z = pd.DataFrame()

            # adding values to final network df
            if not mode_sta_en.empty:
                df_sta2 = pd.DataFrame(
                    data=mode_sta_en.values.reshape(1, len(mode_sta_en)),
                    columns=per,
                    index=[station],
                )
                df_sta = df_sta.append(df_sta2)
            if not dfz_chan.empty:
                dfz_sta2 = pd.DataFrame(
                    data=mode_sta_z.values.reshape(1, len(mode_sta_z)),
                    columns=per,
                    index=[station],
                )
                dfz_sta = dfz_sta.append(dfz_sta2)

        # Add stats, sta, z sta and nlnm to output
        df_dict = {}
        df_dict["stats"] = stats
        df_dict["df_sta"] = df_sta
        df_dict["dfz_sta"] = dfz_sta
        df_dict["nlnm"] = nlnm

    # TODO: for refactor, add check to ensure isinstance is Database
    else:
        database = st
        df_dict = unique_stations_psds_stats(database, network, station, channel)

    # Calculate model values
    out = calc_model_vals(df_dict, network)

    # If plot true, then plot enModel
    if plot:
        plt.close()
        if len(out["enModel"][0]) > 0:
            plt.semilogx(
                out["enModel"][0],
                out["enModel"][1],
                "orangered",
                label="E/N Network Model",
            )

        # Plot z Model
        if len(out["zModel"][0]) > 0:
            plt.semilogx(out["zModel"][0], out["zModel"][1], "c", label="Z Network Model")

        plt.xlabel("Period")
        plt.ylabel("Power(dB)")
        plt.title("Noise Model for Network: " + network)
        plt.grid(True, "both", "both")

        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(list(zip(labels, handles)))
        plt.legend(list(by_label.values()), list(by_label.keys()))

        if ylo is not None and yhi is not None:
            plt.ylim(ylo, yhi)

        if fname is not None:
            name = fname + ".png"
            plt.savefig(fname + ".png", dpi=300)

            update = {"network_noise_image_path": name}
            out.update(update)

        plt.close()

    # If database defined, insert information
    if isinstance(st, Database) or database is not None:
        database.insert_metric(out)
    return out


def _break_large_trace_into_days(tr):
    """Takes in a trace which spans multiple days and creates slices of the trace the size of a single day

    :param tr: Trace of waveform data
    :type tr: obspy.core.trace.Trace
    :return: List of UTCDateTime (start, end) tuples each one day in length
    :rtype: list[tuple(obspy.UTCDateTime, obspy.UTCDateTime)]

    """

    # Initialize output list
    traces_by_day = []
    # Grab start time and end time
    startt = tr.stats.starttime
    endt = tr.stats.endtime
    # Loop through and break out on day boundaries
    while startt.day <= endt.day:
        new_end = UTCDateTime(startt.year, startt.month, startt.day + 1)
        traces_by_day.append((startt, new_end))
        startt = new_end
    return traces_by_day


def _convert_days_to_trace_slices(tr, traces_by_day):
    """
    Takes in a trace wtih a list of date slices and creates a single stream

    :param tr: Trace of waveform data
    :type tr: obspy.core.trace.Trace
    :param traces_by_day: List of UTCDateTime (start, end) tuples each one day in length
    :type trace_by_day: list[tuple(obspy.UTCDateTime, obspy.UTCDateTime)]
    :return: Stream of trace slices, where each trace spans a single day
    :rtype: obspy.core.stream.Stream
    """

    # Initialize output list
    updated_traces = []
    # Loop through traces and create stream object
    for trtups in traces_by_day:
        updated_traces.append(tr.slice(trtups[0], trtups[1], nearest_sample=False))
    return updated_traces


def multiday_traces_to_single_days(st):
    """
    If a stream has a trace which spans multiple days, that trace will be converted
    to multiple traces which span a single day each. Traces which already span less
    than a day are not transformed

    :param st: Stream of waveform data
    :type st: obspy.core.stream.Stream
    """

    # Create empty stream object
    new_st = Stream()
    # Loop through traces in stream and break them apart into traces that span a single day
    for tr in st.traces:
        if (tr.stats.endtime.day - tr.stats.starttime.day) > 1:
            tr_days = _break_large_trace_into_days(tr)
            tr_daily_list = _convert_days_to_trace_slices(tr, tr_days)
            # Add daily trace list to new stream
            new_st += Stream(tr_daily_list)
        else:
            new_st += Stream(tr)
    return new_st


def unique_stations_psds_stats(
    database,
    network,
    station=None,
    channel=None,
    time=None,
    session=None,
    location=None,
    logger=None,
):
    # Initiate empty pandas dataframes
    df_sta = pd.DataFrame()
    dfz_sta = pd.DataFrame()

    # Grab data from psdMetric from databsae
    tb = database.get_metric(
        "psdMetric",
        network=network,
        station=station,
        channel=channel,
        location=location,
        session=session,
        time=time,
    )

    # If returns empty, exit
    if tb.empty:
        logger.error(
            "networkNoiseModel Error: No psdMetric found for specific network in Database, please run psdMetric first"
        )
        return

    # Get list of unique stations
    sta = list(tb.station.unique())
    # Loop through stations
    for i in range(len(sta)):
        station = np.unique(tb.station)[i]
        # Get psds
        tb_psds = tb.uncorrected_psds[tb.station == station]
        chan_psds = []
        for j in range(len(tb_psds)):
            psd = database.extract_masks(tb_psds.iloc[j])
            chan_psds.append(psd)
        stats = psdStatistics(chan_psds, evalresp=None, logger=logger)
        # get period where freq < 10Hz
        try:
            freq = np.asarray(chan_psds[0][0][0])
            per = 1 / freq
        except IndexError:
            # del sta_n[i]
            continue
        # get NLNM where freq < 10Hz
        # get NLNM where freq < 10Hz
        nlnm = stats[0]["nlnm"][np.where(per)]
        # get channels for df  index
        chan = []
        # extract each channel modes
        df_chan = pd.DataFrame()
        dfz_chan = pd.DataFrame()
        for j in range(len(stats)):
            # getting mode for each channel
            try:
                mode = stats[j]["mode"][np.where(per)]
            except IndexError:
                continue
            chan.append(getChannelName(stats[j]["snclq"]))
            # adding modes to df
            df2 = pd.DataFrame(data=mode.reshape(1, len(mode)), columns=per)
            df_chan = df_chan.append(df2, ignore_index=True)
        # adding indexes
        df_chan["snclq"] = chan

        # creating copy of df, so wont get index error as I am dropping Z channels
        df_chan_copy = df_chan
        # loop through to find z channel

        for k in range(len(df_chan.index)):
            name = df_chan.iloc[k].snclq
            if name.endswith("Z"):
                # make z channel have own df
                df = pd.DataFrame(
                    df_chan.iloc[k].values.reshape(1, len(df_chan.iloc[k].values)),
                    columns=list(df_chan.iloc[k].index),
                )
                dfz_chan = dfz_chan.append(df, ignore_index=True)
                # drop z channel from E/N df
                df_chan_copy = df_chan_copy.drop(k)
        # reassigning the df with dropped Z channels to original df
        df_chan = df_chan_copy

        df_chan = df_chan.rename(df_chan["snclq"])
        df_chan = df_chan.drop(["snclq"], axis=1)

        if not df_chan.empty:
            # getting the mean between E/N station modes
            mode_sta_en = df_chan.mean()
        else:
            mode_sta_en = pd.DataFrame()

        if not dfz_chan.empty:
            dfz_chan = dfz_chan.rename(dfz_chan["snclq"])
            dfz_chan = dfz_chan.drop(["snclq"], axis=1)
            mode_sta_z = dfz_chan.mean()
        else:
            mode_sta_z = pd.DataFrame()

        # adding values to final network df
        if mode_sta_en.empty is False:
            df_sta2 = pd.DataFrame(
                data=mode_sta_en.values.reshape(1, len(mode_sta_en)),
                columns=per,
                index=[station],
            )
            df_sta = df_sta.append(df_sta2)
        if dfz_chan.empty is False:
            dfz_sta2 = pd.DataFrame(
                data=mode_sta_z.values.reshape(1, len(mode_sta_z)),
                columns=per,
                index=[station],
            )
            dfz_sta = dfz_sta.append(dfz_sta2)
    df_dict = {}
    df_dict["stats"] = stats
    df_dict["df_sta"] = df_sta
    df_dict["dfz_sta"] = dfz_sta
    df_dict["nlnm"] = nlnm
    return df_dict


def calc_model_vals(df_dict, network):
    """
    Calculate model values
    """

    # Initialize empty lists
    noise_matrix = []
    mode_list = []
    snclq_list = []

    # adding mode, noise matrix, and snclq data to list for out dict
    noise_tmp = []
    mode_tmp = []
    snclq_tmp = []
    # Loop through stats and grab out noise matrix frequency, mode, and snclq
    for j in range(len(df_dict["stats"])):
        noise_tmp.append(df_dict["stats"][j]["noise_matrix_frequency"])
        mode_tmp.append(df_dict["stats"][j]["mode"])
        snclq_tmp.append(df_dict["stats"][j]["snclq"])

    # Append above to initialized lists
    noise_matrix.append(noise_tmp)
    mode_list.append(mode_tmp)
    snclq_list.append(snclq_tmp)

    # TODO: A network name must be provided as a parameter so shape of numpy array is same as pandas dataFrame.
    # If providing network name, networkNoiseModel won't try to perform calculations on different networks
    # As a fix to this, maybe pass in a dict of {'db': Database, 'network': str}

    if df_dict["df_sta"].empty is False:
        # Get difference between station values and nlnm
        df_diff = df_dict["df_sta"] - df_dict["nlnm"]
        # Remove stations with mode value differences >-80 db
        df_diff_bool = df_diff[df_diff < -80].any(axis=1)
        # Convert back to a dataframe
        df_diff_bool = df_diff_bool.to_frame()
        # Get index where values = True
        inds = np.where(df_diff_bool)
        # If inds is not empty, then drop rows in df_sta where df_diff_bool is True using index
        if inds[0]:
            df_dict["df_sta"] = df_dict["df_sta"].drop(inds[0])

    if df_dict["dfz_sta"].empty is False:
        # Get difference between station values and nlnm
        dfz_diff = df_dict["dfz_sta"] - df_dict["nlnm"]
        # Remove stations with mode value differences >-80 db
        dfz_diff_bool = dfz_diff[dfz_diff < -80].any(axis=1)
        # Convert back to a dataframe
        dfz_diff_bool = dfz_diff_bool.to_frame()
        # Get index where values = True
        inds = np.where(dfz_diff_bool)
        # If inds is not empty then drop rows in dfz_sta where dfz_diff_bool is True using the index
        if inds[0]:
            df_dict["dfz_sta"] = df_dict["dfz_sta"].drop(inds[0])

    # Create E/N network noise model ("enModel") and Z network noise model ("zModel") based on mean of remaining
    # stations within output dictionary
    df_sta_out = []
    df_sta_out.append(list(df_dict["df_sta"].columns))
    df_sta_out.append(list(df_dict["df_sta"].mean()))

    dfz_sta_out = []
    dfz_sta_out.append(list(df_dict["dfz_sta"].columns))
    dfz_sta_out.append(list(df_dict["dfz_sta"].mean()))

    # Create output dictionary with E/N and Z network noise models
    out = {
        "enModel": df_sta_out,
        "zModel": dfz_sta_out,
        "stations_en": list(df_dict["df_sta"].index),
        "stations_z": list(df_dict["dfz_sta"].index),
        "metric_name": "networkNoiseModel",
        "network": network,
        "psdStats": {
            "noise_matrix_frequency": noise_matrix,
            "mode": mode_list,
            "snclq": snclq_list,
        },
    }
    return out
