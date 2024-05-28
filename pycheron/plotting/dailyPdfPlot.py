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

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import obspy
import pandas as pd


from pycheron.db.sqllite_db import Database
from pycheron.psd.noise.networkNoiseModel import networkNoiseModel
from pycheron.psd.psdList import psdList
from pycheron.psd.psdStatistics import psdStatistics
from pycheron.util.logger import Logger
from pycheron.util.masks import consecutive
from pycheron.util.format import parse_snclq
from pycheron.util.getMPeriod import get_m_period
from pycheron.metrics.psdMetric import psdMetric


__all__ = [
    "dailyPdfplots",
    "_getDays",
    "_find_nearest",
    "noiseDiff",
    "plot_grid",
    "_plot_colorline",
    "get_pdf_plot_data",
    "plot_grid_data_fill_in",
]


def dailyPdfplots(
    st,
    model="nlnm",
    f_name_grid=None,
    f_name_line=None,
    per_arr=None,
    net_dict=None,
    diff_threshold=50,
    micro_threshold=10,
    banded_threshold=5,
    station=None,
    network=None,
    channel=None,
    session=None,
    location=None,
    database=None,
    logger=None,
):
    """
    Plots a color grid of the daily PDF Mode differences from noise model (nlnm, gsn, or network). Colorgrid plots show
    the magnitude of deviation (in decibels (dB)) of a station’s daily PSD PDF mode from the New Low Noise Model (NLNM)
    (Peterson, 1993), the Global Seismic Network 2004 model (obtained from ported code base), or a user
    provided/calculated network noise model.

    :param st: Obspy stream object or Database Object
    :type st: obspy.core.stream.Stream or pycheron.db.sqllite_db.Database
    :param model: Noise Model to use in differencing from the PDF mode. Options are "gsn", "nlnm", or "network".
                  Using "network" calls networkNoiseModel.py and creates one on the fly or uses input dictionary
                  from the output of a pre-existing network noise model.
    :type model: str
    :param f_name_grid: output file name for png colorgrid plot
    :type f_name_grid: str
    :param f_name_line: output file name for png pdf mode timeline plots
    :type f_name_line: str
    :param per_arr: Period array to use with the models
    :type per_arr: numpy.ndarray
    :param net_dict: Network noise model formatted in same way as output from networkNoiseModel.py
    :type net_dict: dict
    :param diff_threshold: threshold to use to trigger frequency masks.
    :type diff_threshold: int
    :param micro_threshold: Threshold to use to trigger masks when microseism isn't recording healthy energy.
    :type micro_threshold: int
    :param banded_threshold: threshold to use to trigger masks when banded character in colorgrid present.
    :type banded_threshold: int
    :param network: If using database, must supply the network
    :type network: str
    :param station: If using database, must supply the station
    :type station: str
    :param channel: If using database, and want to plot specific channels, if None, all channels
                        associated with network, station, and session will be plotted
    :type channel: str
    :param session: specific session in database
    :type session: str
    :param database: Database object to insert result. If :param st is a Database object, this is not needed.
                     However if it is a Obspy Stream and you would like to save to the DB or read existing data
                     set this parameter to a database object (Default = None)
    :type database: pycheron.db.sqllite_db.Database
    :param logger: Logger object
    :type logger: pycheron.util.logger.Logger

    :return: list of dictionaries with the following keys and types:

        * start_time (`str`): start time of the trace object. Applies to all masks noted below.
        * end_time (`str`): end time of the trace object. Applies to all masks noted below.
        * snclq (`str`): station network channel location quality indicator for the obspy object (e.g., UU.MSU..EHZ)
        * noise_masks (dictionary): dictionary containing the following information:

          * frequency_start (`list`): list of beginning frequency value band where issue(s) start
          * frequency_end (`list`): list of ending frequency value band where issue(s) end

        * microseism_masks (`boolean`): boolean indicating whether the microseism had healthy energy or not
        * banded_masks (`boolean`): boolean indicated whether banded issue was present, True indicates likely unhealthy
          station

    :rtype: list

    Colorgrid plots, daily PDF mode timeline plots, and station ranking plots were first established for use in the
    USArray TA project (Busby et. al., 2018) and code was ported from Perl and Hypertext Preprocessor (PHP) to Python.
    Code was adapted and augmented for use within Pycheron.

    #. Busby, R. W., R.L. Woodward, K.A. Hafner, F. L. Vernon, and A.M. Frasseto (2018). The Design and Implementation
    of EarthScope's USArray Transportable Array in the Conterminous United States and Southern Canada (Rep.).

    **Parameter Notes**

    * diff_threshold (`int`):
          * Threshold to use to trigger frequency masks. This is the difference from the nlnm,
            gsn, or network noise models desired to trigger a mask. DEFAULT = 50 and is probably
            sufficient for most cases as this will be the red color and above within the color
            grid plots. This is used as a greater than or equal to threshold.
    * micro_threshold (`int`):
          * Threshold to use to trigger masks when microseism isn't recording healthy energy.
            This threshold is based on difference from nlnm, gsn, or network noise models desired
            to trigger a mask. DEFAULT = 10, so the two darkest blues on the colorgrid plots.
            This will be utilized as a less than or equal to threshold.
    * banded_threshold (`int`):
          * threshold to use to trigger masks when banded character in colorgrid present. This
            often indicates an issue with the station, often unhealthy behavior, whether a dead
            channel, noise issues, or other equipment issues. This trigger is to let the user
            know there's an issue with the health of the station. This threshold is based on
            the average difference across each element in the array (or rows of the dataFrame)
            Lower numbers indicate there isn't much difference between the elements, which
            suggests an unhealthy station, while larger numbers suggest a healthier change
            as would be expected in a properly functioning station. DEFAULT = 5, and is probably
            a good threshold for the majority of cases.


    **Example**

    .. code-block:: python

       # Import libs
       from obspy.clients.fdsn import Client
       from pycheron.plotting.dailyPdfPlot import dailyPdfplots
       from obspy import UTCDateTime
       client = Client("IRIS")

       # Define start/end time and then get stream object
       starttime = UTCDateTime("2012-12-12T00:00:00.000")
       endtime = UTCDateTime("2012-12-16T00:00:00.000")
       st = client.get_waveforms("AK","GHO","","BHE",starttime, endtime)
       dailyPdfplots(st, model="gsn")

    .. image:: _static/dailyPDFPlotGrid.png

    .. image:: _static/dailyPDFPlotLine.png
    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # if st is a path to directory. open, read in psds and add to list
    if isinstance(st, Database):
        data = _pdf_plots_from_database(
            st,
            model=model,
            f_name_grid=f_name_grid,
            f_name_line=f_name_line,
            per_arr=per_arr,
            net_dict=None,
            diff_threshold=50,
            micro_threshold=10,
            banded_threshold=5,
            station=station,
            network=network,
            channel=channel,
            session=session,
            location=location,
            logger=logger,
        )
        return data

    # Otherwise loop through stream object
    else:
        newStreams = []
        # loop through traces in stream
        for i in range(len(st)):
            tr = st[i]
            # if trace is longer than 24 hours, split traces, create new stream.
            # newStream is a list containing 1 stream per channel. Each stream holds 1 trace per day
            if tr.stats.endtime - tr.stats.starttime >= 86400:  # secs in a day
                stN = _getDays(tr)
                newStreams.append(stN)
            else:
                stN = obspy.Stream(traces=[tr])
                newStreams.append(stN)
        # read in gsn_model. GSN model obtained from ported code base
        if model == "gsn":
            # Create empty arrays to fill for thresholds
            data = []
            period_masks = []
            freq_start = []
            freq_end = []
            # Grab out gsn noise model file and read into dataFrame
            gsn = os.path.dirname(__file__) + "/gsn2004pdf_interpolated.txt"
            model_gsn = pd.read_table(gsn, names=["period", "z_power", "h_power"], delimiter=r"\s+")
            # Get period and power from model
            m_period = model_gsn["period"]
            m_power = model_gsn["h_power"]

            # Loop through traces in stream and parse out snclq information
            for stream in newStreams:
                tr = stream[0]
                network, station, channel, location, quality = parse_snclq(stream[0].get_id())

                # If database object, check to see if we already have the psds calculated
                if database is not None:
                    time = (
                        tr.stats.starttime.isoformat(),
                        tr.stats.endtime.isoformat(),
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

                    # if returns no data, calculate and insert into db
                    # Include check to ensure length of values from the db is same as number of streams entered
                    if tb.empty or len(tb.index) != len(stream):
                        psds = psdList(stream)

                    else:
                        tb_psds = tb.uncorrected_psds[tb.channel == channel]
                        psds = []
                        for j in range(len(tb_psds)):
                            psd = database.extract_masks(tb_psds.iloc[j])
                            psds.append(psd)

                # If not in database, calculate the psds
                else:
                    psds = psdList(stream)
                # Calculate statistics
                stats = psdStatistics(psds, logger=logger, database=database)

                # If Z channel, compare to right power
                if channel.endswith("Z"):
                    m_power = model_gsn["z_power"]

                # Create dataframes to fill
                df_line = pd.DataFrame()
                df = pd.DataFrame()
                date_index = []
                # loop through traces in NEW stream
                for j in range(len(stream)):
                    pmode = stats[j]["mode"]
                    date = str(stream[j].stats.starttime)
                    freq = np.asarray(psds[j][0][0])
                    period = 1 / freq

                    date_index.append(date[0:10])
                    # loop through all periods and find nearest in gsn model
                    diff, per = noiseDiff(period, pmode, m_power, m_period)

                    df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                    df = pd.concat([df, df2], ignore_index=True)

                    df_line2 = pd.DataFrame(data=pmode.reshape(1, len(pmode)), columns=per)
                    df_line = pd.concat([df_line, df_line2], ignore_index=True)

                # Drop and rename
                df["Date"] = date_index
                df = df.rename(df["Date"])
                df = df.drop(["Date"], axis=1)

                df_line["Date"] = date_index
                df_line = df_line.rename(df_line["Date"])
                df_line = df_line.drop(["Date"], axis=1)

                # Drop digitizer values
                df = df.drop([df.columns[-1]], axis=1)
                df = df.drop([df.columns[-2]], axis=1)
                df = df.drop([df.columns[-3]], axis=1)

                snclq = stream[0].get_id()

                # Plot colorgrid and mode timelines
                update_grid = plot_grid(df, snclq, plot="dailyPDF", f_name=f_name_grid)
                update_line = _plot_colorline(df_line, snclq, per_arr, f_name=f_name_line)

                # dropping digitizer values
                try:
                    ### Trigger based on large differences from noise model ###
                    # Trigger mask output based on diff_threshold value -- looks for
                    thresh = np.where(df >= diff_threshold)[1]
                    # Get consecutive indices out

                    thresh_cons = consecutive(thresh, stepsize=1)

                    # Loop through consecutive indices and find where freq values where threshold exceeds those values
                    for k in range(len(thresh_cons)):
                        period_masks.append(df.columns[thresh_cons[k]])

                    # Grab out period masks and aggregate them into start/stop frequencies
                    for k in range(len(period_masks)):
                        freq_start.append(1 / float(period_masks[k][0]))
                        freq_end.append(1 / float(period_masks[k][-1]))

                    ### Trigger based on unhealthy differences from noise model in the microseism band ###
                    # Microseism band is 0.1 - 1.0 Hz, which translates to 1 - 10s periods
                    # Find indices that are between 1-10s period so that we can subset the dataframe accordingly
                    # First convert columns to list, then convert strings to floats
                    cols = list(df.columns)
                    cols_f = [float(i) for i in cols]
                    micro_dead = False
                    # Find index of max period (10s) and min period (1s)
                    max_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 10))
                    min_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 1))

                    # Find where values are below threshold, indicating a potentially dead channel
                    micro_thresh = np.where(df.iloc[0][max_period:min_period] <= micro_threshold)

                    # Check if length of micro_thresh equal to number of elements in period range of microseism, this
                    # means that all of them are less than or equal to the threshold, indicating a dead channel
                    if len(micro_thresh[0]) == len(df.iloc[0][max_period:min_period]):
                        micro_dead = True

                    # Also check if the mean is less than or equal to micro_threshold, indicating a dead channel too,
                    # but perhaps not all elements are less than the micro_threshold
                    if np.mean(df.iloc[0][max_period:min_period]) <= micro_threshold:
                        micro_dead = True

                    ### Trigger on banded character of colorgrid indicating issue ###
                    # This could be a dead channel, noise issues, or other issues.
                    # This trigger is to let the user know that there's an issue with the health of the station.
                    # When we add more complexity to these metrics and review more colorgrids,
                    # we might be able to get more specific here ###

                    # Get rows of dataFrame into a list
                    rows = list(df.iloc[0])
                    banded_issue = False

                    # Then find the average of the difference. Low averages generally indicate that a banded character
                    # is present in the colorGrid plot, while high averages are more likely to indicate that the station
                    # is healthy
                    banded = np.mean(np.diff(rows))

                    if banded <= banded_threshold:
                        banded_issue = True

                    # Create dictionary for every stream object, then append to data to encompass all information per
                    # stream
                    d = {
                        "start_time": tr.stats.starttime.isoformat(),
                        "end_time": tr.stats.starttime.isoformat(),
                        "snclq": tr.get_id(),
                        "noise_masks": {
                            "frequency_start": freq_start,
                            "frequency_end": freq_end,
                        },
                        "microseism_masks": micro_dead,
                        "banded_masks": banded_issue,
                        "metric_name": "dailyPdfPlot",
                    }
                    d.update(update_grid)
                    d.update(update_line)
                    data.append(d)
                # If error out, then set up output dictionary object
                except IndexError:
                    d = {
                        "metric_name": "dailyPdfPlot",
                        "start_time": tr.stats.starttime,
                        "end_time": tr.stats.endtime,
                        "snclq": snclq,
                    }
                    d.update(update_grid)  # add in the plot image and keys
                    d.update(update_line)
                    data.append(d)
                    print("dailyPDFPlot Info: No threshold created.")
                    continue

        # If model nlnm
        elif model == "nlnm":
            # Initialize empty arrays for filling in periods/frequencies and also to fill dictionary
            data = []
            period_masks = []
            freq_start = []
            freq_end = []
            # Loop through traces and parse out snclq object
            for stream in newStreams:
                tr = stream[0]

                network, station, channel, location, quality = parse_snclq(stream[0].get_id())

                # If database object, check to see if we already have the psds calculated
                if database is not None:
                    time = (
                        tr.stats.starttime.isoformat(),
                        tr.stats.endtime.isoformat(),
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

                    # if returns no data, calculate and insert into db
                    if tb.empty or len(tb.index) != len(stream):
                        psds = psdList(stream)
                        psds_metrics_for_db = psdMetric(stream)
                        database.insert_metric(psds_metrics_for_db)
                        freq = psds[0][0][0]
                        period = 1 / freq

                    else:
                        tb_psds = tb.uncorrected_psds[tb.channel == channel]
                        psds = []
                        for j in range(len(tb_psds)):
                            psd = database.extract_masks(tb_psds.iloc[j])
                            psds.append(psd)

                        freq = np.asarray(psds[0][0][0])
                        period = 1 / freq
                # If database not defined, then calculate psds, grab out freq and period
                else:
                    psds = psdList(stream)
                    freq = psds[0][0][0]
                    period = 1 / freq
                # Calculate statistics
                stats = psdStatistics(psds, logger=logger, database=database)
                # Create empty dataframes to fill
                df_line = pd.DataFrame()
                df = pd.DataFrame()
                date_index = []
                # loop through traces in NEW stream
                for val in range(len(stream)):
                    pmode = stats[val]["mode"]
                    m_power = stats[val]["nlnm"]
                    # m_period = 1 / stats[j]["noise_matrix_frequency"]
                    # TODO: Should this be a np array? or a list? Originally a list, but functionality was broken
                    m_period = [1.0 / x for x in stats[val]["noise_matrix_frequency"]]
                    date = str(stream[val].stats.starttime)
                    date_index.append(date[0:10])
                    # loop through all periods and find nearest in nlnm model
                    diff, per = noiseDiff(period, pmode, m_power, m_period)

                    df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                    df = pd.concat([df, df2], ignore_index=True)

                    df_line2 = pd.DataFrame(data=pmode.reshape(1, len(pmode)), columns=per)
                    df_line = pd.concat([df_line, df_line2], ignore_index=True)

                # Drop and rename
                df["Date"] = date_index
                df = df.rename(df["Date"])
                df = df.drop(["Date"], axis=1)

                df_line["Date"] = date_index
                df_line = df_line.rename(df_line["Date"])
                df_line = df_line.drop(["Date"], axis=1)

                # Drop digitizer values
                df = df.drop([df.columns[-1]], axis=1)
                df = df.drop([df.columns[-2]], axis=1)

                snclq = stream[0].get_id()

                # Plot colorgrid and mode timelines
                update_grid = plot_grid(df, snclq, plot="dailyPDF", f_name=f_name_grid)
                update_line = _plot_colorline(df_line, snclq, per_arr, f_name=f_name_line)
                try:
                    ### Trigger based on large differences from noise model ###
                    # Trigger mask output based on diff_threshold value -- looks for
                    thresh = np.where(df >= diff_threshold)[1]
                    # Get consecutive indices out
                    thresh_cons = consecutive(thresh, stepsize=1)

                    # Loop through consecutive indices and find where freq values where threshold exceeds those values
                    for k in range(len(thresh_cons)):
                        period_masks.append(df.columns[thresh_cons[k]])

                    # Grab out period masks and aggregate them into start/stop frequencies
                    for k in range(len(period_masks)):
                        freq_start.append(1 / float(period_masks[k][0]))
                        freq_end.append(1 / float(period_masks[k][-1]))

                    ### Trigger based on unhealthy differences from noise model in the microseism band ###
                    # Microseism band is 0.1 - 1.0 Hz, which translates to 1 - 10s periods
                    # Find indices that are between 1-10s period so that we can subset the dataframe accordingly
                    # First convert columns to list, then convert strings to floats
                    cols = list(df.columns)
                    cols_f = [float(i) for i in cols]
                    micro_dead = False
                    # Find index of max period (10s) and min period (1s)
                    max_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 10))
                    min_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 1))

                    # Find where values are below threshold, indicating a potentially dead channel
                    micro_thresh = np.where(df.iloc[0][max_period:min_period] <= micro_threshold)

                    # Check if length of micro_thresh equal to number of elements in period range of microseism, this
                    # means that all of them are less than or equal to the threshold, indicating a dead channel
                    if len(micro_thresh[0]) == len(df.iloc[0][max_period:min_period]):
                        micro_dead = True

                    # Also check if the mean is less than or equal to micro_threshold, indicating a dead channel too,
                    # but perhaps not all elements are less than the micro_threshold
                    if np.mean(df.iloc[0][max_period:min_period]) <= micro_threshold:
                        micro_dead = True

                    ### Trigger on banded character of colorgrid indicating issue ###
                    # This could be a dead channel, noise issues, or other issues.
                    # This trigger is to let the user know that there's an issue with the health of the station.
                    # When we add more complexity to these metrics and
                    # review more colorgrids we might be able to get more specific here ###

                    # Get rows of dataFrame into a list
                    rows = list(df.iloc[0])
                    banded_issue = False

                    # Then find the average of the difference.
                    # Low averages generally indicate that a banded character is
                    # present in the colorGrid plot, while high averages are more likely to indicate that
                    # the station is healthy
                    banded = np.mean(np.diff(rows))

                    if banded <= banded_threshold:
                        banded_issue = True

                    # Create dictionary for every stream object,
                    # then append to data to encompass all information per stream
                    d = {
                        "start_time": tr.stats.starttime.isoformat(),
                        "end_time": tr.stats.starttime.isoformat(),
                        "snclq": tr.get_id(),
                        "noise_masks": {
                            "frequency_start": freq_start,
                            "frequency_end": freq_end,
                        },
                        "microseism_masks": micro_dead,
                        "banded_masks": banded_issue,
                        "metric_name": "dailyPdfPlot",
                    }
                    d.update(update_grid)
                    d.update(update_line)
                    data.append(d)
                # If error out, then set up output dictionary object
                except IndexError:
                    d = {
                        "metric_name": "dailyPdfPlot",
                        "start_time": tr.stats.starttime,
                        "end_time": tr.stats.endtime,
                        "snclq": snclq,
                    }
                    d.update(update_grid)  # add in the plot image and keys
                    d.update(update_line)
                    data.append(d)
                    print("dailyPDFPlot Info: No threshold created.")
                    continue

        # read in network model if available
        elif model == "network":
            # Create empty arrays to fill for thresholds
            data = []
            period_masks = []
            freq_start = []
            freq_end = []
            # Create variables for period, enModel, zModel
            if net_dict is not None:
                # Grab out z_power, en_power models and period
                z_power = net_dict["zModel"][1]
                en_power = net_dict["enModel"][1]
                # Get indices of period from the zModel, then convert to a series
                m_period = get_m_period(net_dict)
                #
            else:
                # If a previous dictionary is not defined, create a network noise model on the fly. Recommended that
                # we have many stations in the stream object and > 1 week of data for this
                net_noise = networkNoiseModel(st)
                # Grab out z_power, en_power models and period
                z_power = net_noise["zModel"][1]
                en_power = net_noise["enModel"][1]
                # Get indices of period from the zModel, then convert to a series
                m_period = get_m_period(net_noise)

            # Loop through traces and parse out snclq information
            for stream in newStreams:
                tr = stream[0]
                network, station, channel, location, quality = parse_snclq(stream[0].get_id())

                # If database object, check to see if we already have the psds calculated
                if database is not None:
                    time = (
                        tr.stats.starttime.isoformat(),
                        tr.stats.endtime.isoformat(),
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

                    # if returns no data, calculate and insert into db
                    if tb.empty or len(tb.index) != len(stream):
                        psds = psdList(stream)
                        psds_metrics_for_db = psdMetric(stream)
                        database.insert_metric(psds_metrics_for_db)
                        freq = psds[0][0][0]
                        period = 1 / freq

                    else:
                        tb_psds = tb.uncorrected_psds[tb.channel == channel]
                        psds = []
                        for j in range(len(tb_psds)):
                            psd = database.extract_masks(tb_psds.iloc[j])
                            psds.append(psd)
                # If database not defined, calculate psds
                else:
                    psds = psdList(stream)
                # Calculate psdStatistics
                stats = psdStatistics(psds, logger=logger, database=database)

                # Determine the channel, so we determine which network noise model to use, Z or EN
                if channel.endswith("Z"):
                    m_power = z_power
                else:
                    m_power = en_power
                # Set up empty dataframes to fill
                df_line = pd.DataFrame()
                df = pd.DataFrame()
                date_index = []
                # loop through traces in NEW stream
                for j in range(len(stream)):
                    pmode = stats[j]["mode"]
                    date = str(stream[j].stats.starttime)
                    freq = np.asarray(psds[j][0][0])
                    period = 1 / freq

                    date_index.append(date[0:10])
                    # loop through all periods and find nearest in network model
                    diff, per = noiseDiff(period, pmode, m_power, m_period)

                    df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                    df = pd.concat([df, df2], ignore_index=True)

                    df_line2 = pd.DataFrame(data=pmode.reshape(1, len(pmode)), columns=per)
                    df_line = pd.concat([df_line, df_line2], ignore_index=True)

                # Drop and rename
                df["Date"] = date_index
                df = df.rename(df["Date"])
                df = df.drop(["Date"], axis=1)

                df_line["Date"] = date_index
                df_line = df_line.rename(df_line["Date"])
                df_line = df_line.drop(["Date"], axis=1)
                # Drop digitizer values
                df = df.drop([df.columns[-1]], axis=1)
                df = df.drop([df.columns[-2]], axis=1)
                df = df.drop([df.columns[-3]], axis=1)

                snclq = stream[0].get_id()

                # Plot colorgrid and mode timelines
                update_grid = plot_grid(df, snclq, plot="dailyPDF", f_name=f_name_grid)
                update_line = _plot_colorline(df_line, snclq, per_arr, f_name=f_name_line)

                # dropping digitizer values
                try:
                    ### Trigger based on large differences from noise model ###
                    # Trigger mask output based on diff_threshold value -- looks for
                    thresh = np.where(df >= diff_threshold)[1]
                    # Get consecutive indices out
                    thresh_cons = consecutive(thresh, stepsize=1)

                    # Loop through consecutive indices and find where freq values where threshold exceeds those values
                    for k in range(len(thresh_cons)):
                        period_masks.append(df.columns[thresh_cons[k]])

                    # Grab out period masks and aggregate them into start/stop frequencies
                    for k in range(len(period_masks)):
                        freq_start.append(1 / float(period_masks[k][0]))
                        freq_end.append(1 / float(period_masks[k][-1]))

                    ### Trigger based on unhealthy differences from noise model in the microseism band ###
                    # Microseism band is 0.1 - 1.0 Hz, which translates to 1 - 10s periods
                    # Find indices that are between 1-10s period so that we can subset the dataframe accordingly
                    # First convert columns to list, then convert strings to floats
                    cols = list(df.columns)
                    cols_f = [float(i) for i in cols]
                    micro_dead = False
                    # Find index of max period (10s) and min period (1s)
                    max_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 10))
                    min_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 1))

                    # Find where values are below threshold, indicating a potentially dead channel
                    micro_thresh = np.where(df.iloc[0][max_period:min_period] <= micro_threshold)

                    # Check if length of micro_thresh equal to number of elements in period range of microseism, this
                    # means that all of them are less than or equal to the threshold, indicating a dead channel
                    if len(micro_thresh[0]) == len(df.iloc[0][max_period:min_period]):
                        micro_dead = True

                    # Also check if the mean is less than or equal to micro_threshold, indicating a dead channel too,
                    # but perhaps not all elements are less than the micro_threshold
                    if np.mean(df.iloc[0][max_period:min_period]) <= micro_threshold:
                        micro_dead = True

                    ### Trigger on banded character of colorgrid indicating issue ###
                    # This could be a dead channel, noise issues, or other issues.
                    # This trigger is to let the user know thatthere's an issue with the health of the station.
                    # When we add more complexity to these metrics and
                    # review more colorgrids we might be able to get more specific here ###

                    # Get rows of dataFrame into a list
                    rows = list(df.iloc[0])
                    banded_issue = False

                    # Then find the average of the difference. Low averages generally indicate that a banded character
                    # is present in the colorGrid plot, while high averages are more likely to indicate that the station
                    # is healthy
                    banded = np.mean(np.diff(rows))

                    if banded <= banded_threshold:
                        banded_issue = True

                    # Create dictionary for every stream object, then append to data to encompass all information per
                    # stream
                    d = {
                        "start_time": tr.stats.starttime,
                        "end_time": tr.stats.starttime,
                        "snclq": tr.get_id(),
                        "noise_masks": {
                            "frequency_start": freq_start,
                            "frequency_end": freq_end,
                        },
                        "microseism_masks": micro_dead,
                        "banded_masks": banded_issue,
                        "metric_name": "dailyPdfPlot",
                    }
                    d.update(update_grid)
                    d.update(update_line)
                    data.append(d)
                # If error out, create output dictionary
                except IndexError:
                    d = {
                        "metric_name": "dailyPdfPlot",
                        "start_time": tr.stats.starttime,
                        "end_time": tr.stats.endtime,
                        "snclq": snclq,
                    }
                    d.update(update_grid)  # add in the plot image and keys
                    d.update(update_line)
                    data.append(d)
                    print("dailyPDFPlot Info: No threshold created.")
                    continue
        # Ensure user provides proper model input
        else:
            logger.error(
                "dailyPDFplots(): Error model: {m} is not a valid model, please choose from nlnm, nnm, or gsn".format(
                    m=model
                )
            )

        # Insert data into database if available
        if database is not None:
            database.insert_metric(data)

        return data


def get_pdf_plot_data(
    database,
    network,
    station,
    location=None,
    channel=None,
    session=None,
    model="nlnm",
    f_name_grid=None,
    f_name_line=None,
    per_arr=None,
    net_dict=None,
    diff_threshold=50,
    micro_threshold=10,
    banded_threshold=5,
    logger=None,
):

    if model == "nlnm":

        tb = database.get_metric(
            "psdMetric",
            network=network,
            station=station,
            channel=channel,
            location=location,
            session=session,
        )

        if tb.empty:
            logger.error(
                "dailyPDFPlot: Error no PSDS found in DB for Network: \
                    {n}, Station: {s}, Channel: {c}, Session: {se}".format(
                    n=network, s=station, c=channel, se=session
                )
            )
            return

        for i in range(len(tb.channel.unique())):
            df_line = pd.DataFrame()
            df = pd.DataFrame()
            date_index = []

            channel = tb.channel.iloc[i]
            tb_psds = tb.uncorrected_psds[tb.channel == channel]
            chan_psds = []
            for j in range(len(tb_psds)):
                psd = database.extract_masks(tb_psds.iloc[j])
                chan_psds.append(psd)
            stats = psdStatistics(chan_psds, logger=logger, database=database)

            df = pd.DataFrame()
            df_line = pd.DataFrame()
            # loop through traces in NEW stream
            for j in range(len(chan_psds)):
                pmode = stats[j]["mode"]
                m_power = np.asarray(stats[j]["nlnm"])
                m_period = 1 / np.asarray(stats[j]["noise_matrix_frequency"])
                snclq = tb.snclq[tb.channel == channel].iloc[j]
                date = tb.start_time[tb.channel == channel].iloc[j]
                # starttime = tb.start_time[tb.channel == channel].iloc[j]
                # endtime = tb.end_time[tb.channel == channel].iloc[j]
                freq = np.asarray(chan_psds[j][0][0])
                period = 1 / freq

                date_index.append(date[0:10])
                # loop through all periods and find nearest in gsn model
                diff, per = noiseDiff(period, pmode, m_power, m_period)

                df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                df = pd.concat([df, df2], ignore_index=True)

                df_line2 = pd.DataFrame(data=pmode.reshape(1, len(pmode)), columns=period)
                df_line = pd.concat([df_line, df_line2], ignore_index=True)

            df["Date"] = date_index
            df = df.rename(df["Date"])
            df = df.drop(["Date"], axis=1)

            df_line["Date"] = date_index
            df_line = df_line.rename(df_line["Date"])
            df_line = df_line.drop(["Date"], axis=1)

            # dropping digitizer values
            df = df.drop([df.columns[-1]], axis=1)
            df = df.drop([df.columns[-2]], axis=1)

            return df, df_line, snclq


def _pdf_plots_from_database(
    database,
    network,
    station,
    location=None,
    channel=None,
    session=None,
    model="nlnm",
    f_name_grid=None,
    f_name_line=None,
    per_arr=None,
    net_dict=None,
    diff_threshold=50,
    micro_threshold=10,
    banded_threshold=5,
    logger=None,
):
    #import pdb; pdb.set_trace()
    if model == "gsn":
        # Create empty arrays to fill for thresholds
        data = []
        period_masks = []
        freq_start = []
        freq_end = []
        # Grab out gsn noise model file and read into dataFrame
        gsn = os.path.dirname(__file__) + "/gsn2004pdf_interpolated.txt"
        model_gsn = pd.read_table(gsn, names=["period", "z_power", "h_power"], delimiter=r"\s+")
        m_period = model_gsn["period"]
        m_power = model_gsn["h_power"]

        tb = database.get_metric(
            "psdMetric",
            network=network,
            station=station,
            channel=channel,
            location=location,
            session=session,
        )

        if tb.empty:
            logger.error(
                "dailyPDFPlot: Error no PSDS found in DB for Network: \
                    {n}, Station: {s}, Channel: {c}, Session: {se}".format(
                    n=network, s=station, c=channel, se=session
                )
            )
            return

        for i in range(len(tb.channel.unique())):
            df_line = pd.DataFrame()
            df = pd.DataFrame()
            date_index = []

            channel = tb.channel.iloc[i]
            tb_psds = tb.uncorrected_psds[tb.channel == channel]
            chan_psds = []
            for j in range(len(tb_psds)):
                psd = database.extract_masks(tb_psds.iloc[j])
                chan_psds.append(psd)
            stats = psdStatistics(chan_psds, logger=logger, database=database)

            if channel.endswith("Z"):
                m_power = model_gsn["z_power"]

            for k in range(len(chan_psds)):
                pmode = stats[k]["mode"]
                snclq = tb.snclq[tb.channel == channel].iloc[k]
                date = tb.start_time[tb.channel == channel].iloc[k]
                starttime = tb.start_time[tb.channel == channel].iloc[k]
                endtime = tb.end_time[tb.channel == channel].iloc[k]
                freq = np.asarray(chan_psds[k][0][0])
                period = 1 / freq

                date_index.append(date[0:10])
                # loop through all periods and find nearest in gsn model
                diff, per = noiseDiff(period, pmode, m_power, m_period)

                df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                df = pd.concat([df, df2], ignore_index=True)

                df_line2 = pd.DataFrame(data=pmode.reshape(1, len(pmode)), columns=per)
                df_line = pd.concat([df_line, df_line2], ignore_index=True)

            df["Date"] = date_index
            df = df.rename(df["Date"])
            df = df.drop(["Date"], axis=1)

            df_line["Date"] = date_index
            df_line = df_line.rename(df_line["Date"])
            df_line = df_line.drop(["Date"], axis=1)

            # dropping digitizer values
            print(f"Network: {network}, Station: {station}, Channel: {channel}")
            df = df.drop([df.columns[-1]], axis=1)
            df = df.drop([df.columns[-2]], axis=1)
            df = df.drop([df.columns[-3]], axis=1)

            update_grid = plot_grid(df, snclq, plot="dailyPDF", f_name=f_name_grid)
            update_line = _plot_colorline(df_line, snclq, per_arr, f_name=f_name_line)

            try:
                ### Trigger based on large differences from noise model ###
                # Trigger mask output based on diff_threshold value -- looks for
                thresh = np.where(df >= diff_threshold)[1]
                # Get consecutive indices out
                thresh_cons = consecutive(thresh, stepsize=1)

                # Loop through consecutive indices and find where freq values where threshold exceeds those values
                for k in range(len(thresh_cons)):
                    period_masks.append(df.columns[thresh_cons[k]])

                # Grab out period masks and aggregate them into start/stop frequencies
                for k in range(len(period_masks)):
                    freq_start.append(1 / float(period_masks[k][0]))
                    freq_end.append(1 / float(period_masks[k][-1]))

                ### Trigger based on unhealthy differences from noise model in the microseism band ###
                # Microseism band is 0.1 - 1.0 Hz, which translates to 1 - 10s periods
                # Find indices that are between 1-10s period so that we can subset the dataframe accordingly
                # First convert columns to list, then convert strings to floats
                cols = list(df.columns)
                cols_f = [float(i) for i in cols]
                micro_dead = False
                # Find index of max period (10s) and min period (1s)
                max_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 10))
                min_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 1))

                # Find where values are below threshold, indicating a potentially dead channel
                micro_thresh = np.where(df.iloc[0][max_period:min_period] <= micro_threshold)

                # Check if length of micro_thresh equal to number of elements in period range of microseism, this
                # means that all of them are less than or equal to the threshold, indicating a dead channel
                if len(micro_thresh[0]) == len(df.iloc[0][max_period:min_period]):
                    micro_dead = True

                # Also check if the mean is less than or equal to micro_threshold, indicating a dead channel too,
                # but perhaps not all elements are less than the micro_threshold
                if np.mean(df.iloc[0][max_period:min_period]) <= micro_threshold:
                    micro_dead = True

                ### Trigger on banded character of colorgrid indicating issue ###
                # This could be a dead channel, noise issues, or other issues. This trigger is to let the user know that
                # there's an issue with the health of the station. When we add more complexity to these metrics and
                # review more colorgrids we might be able to get more specific here ###

                # Get rows of dataFrame into a list
                rows = list(df.iloc[0])
                banded_issue = False

                # Then find the average of the difference. Low averages generally indicate that a banded character
                # is present in the colorGrid plot, while high averages are more likely to indicate that the station
                # is healthy
                banded = np.mean(np.diff(rows))

                if banded <= banded_threshold:
                    banded_issue = True

                # Create dictionary for every stream object, then append to data to encompass all information per
                # stream
                d = {
                    "start_time": starttime,
                    "end_time": endtime,
                    "snclq": snclq,
                    "noise_masks": {
                        "frequency_start": freq_start,
                        "frequency_end": freq_end,
                    },
                    "microseism_masks": micro_dead,
                    "banded_masks": banded_issue,
                    "metric_name": "dailyPdfPlot",
                }
                d.update(update_grid)  # add in the plot image and keys
                d.update(update_line)
                data.append(d)

            except (IndexError, ZeroDivisionError):
                d = {
                    "metric_name": "dailyPdfPlot",
                    "start_time": starttime,
                    "end_time": endtime,
                    "snclq": snclq,
                }
                d.update(update_grid)  # add in the plot image and keys
                d.update(update_line)
                data.append(d)
                print(("dailyPDFPlot Info: No threshold created for SNCLQ: " + snclq))
                continue

    if model == "nlnm":
        # Initialize empty arrays for filling in periods/frequencies and also to fill dictionary
        data = []
        period_masks = []
        freq_start = []
        freq_end = []

        tb = database.get_metric(
            "psdMetric",
            network=network,
            station=station,
            channel=channel,
            location=location,
            session=session,
        )

        if tb.empty:
            logger.error(
                "dailyPDFPlot: Error no PSDS found in DB for Network: \
                    {n}, Station: {s}, Channel: {c}, Session: {se}".format(
                    n=network, s=station, c=channel, se=session
                )
            )
            return

        for i in range(len(tb.channel.unique())):
            df_line = pd.DataFrame()
            df = pd.DataFrame()
            date_index = []

            channel = tb.channel.iloc[i]
            tb_psds = tb.uncorrected_psds[tb.channel == channel]
            chan_psds = []
            for j in range(len(tb_psds)):
                psd = database.extract_masks(tb_psds.iloc[j])
                chan_psds.append(psd)
            stats = psdStatistics(chan_psds, logger=logger, database=database)

            df = pd.DataFrame()
            df_line = pd.DataFrame()
            # loop through traces in NEW stream
            for j in range(len(chan_psds)):
                pmode = stats[j]["mode"]
                m_power = np.asarray(stats[j]["nlnm"])
                m_period = 1 / np.asarray(stats[j]["noise_matrix_frequency"])
                snclq = tb.snclq[tb.channel == channel].iloc[j]
                date = tb.start_time[tb.channel == channel].iloc[j]
                starttime = tb.start_time[tb.channel == channel].iloc[j]
                endtime = tb.end_time[tb.channel == channel].iloc[j]
                freq = np.asarray(chan_psds[j][0][0])
                period = 1 / freq

                date_index.append(date[0:10])
                # loop through all periods and find nearest in gsn model
                diff, per = noiseDiff(period, pmode, m_power, m_period)

                df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                df = pd.concat([df, df2], ignore_index=True)

                df_line2 = pd.DataFrame(data=pmode.reshape(1, len(pmode)), columns=period)
                df_line = pd.concat([df_line, df_line2], ignore_index=True)

            df["Date"] = date_index
            df = df.rename(df["Date"])
            df = df.drop(["Date"], axis=1)

            df_line["Date"] = date_index
            df_line = df_line.rename(df_line["Date"])
            df_line = df_line.drop(["Date"], axis=1)

            # dropping digitizer values
            df = df.drop([df.columns[-1]], axis=1)
            df = df.drop([df.columns[-2]], axis=1)

            update_grid = plot_grid(df, snclq, plot="dailyPDF", f_name=f_name_grid)
            update_line = _plot_colorline(df_line, snclq, per_arr, f_name=f_name_line)
            try:
                ### Trigger based on large differences from noise model ###
                # Trigger mask output based on diff_threshold value -- looks for
                thresh = np.where(df >= diff_threshold)[1]
                # Get consecutive indices out
                thresh_cons = consecutive(thresh, stepsize=1)

                # Loop through consecutive indices and find where freq values where threshold exceeds those values
                for k in range(len(thresh_cons)):
                    period_masks.append(df.columns[thresh_cons[k]])

                # Grab out period masks and aggregate them into start/stop frequencies
                for k in range(len(period_masks)):
                    freq_start.append(1 / float(period_masks[k][0]))
                    freq_end.append(1 / float(period_masks[k][-1]))

                ### Trigger based on unhealthy differences from noise model in the microseism band ###
                # Microseism band is 0.1 - 1.0 Hz, which translates to 1 - 10s periods
                # Find indices that are between 1-10s period so that we can subset the dataframe accordingly
                # First convert columns to list, then convert strings to floats
                cols = list(df.columns)
                cols_f = [float(i) for i in cols]
                micro_dead = False
                # Find index of max period (10s) and min period (1s)
                max_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 10))
                min_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 1))

                # Find where values are below threshold, indicating a potentially dead channel
                micro_thresh = np.where(df.iloc[0][max_period:min_period] <= micro_threshold)

                # Check if length of micro_thresh equal to number of elements in period range of microseism, this
                # means that all of them are less than or equal to the threshold, indicating a dead channel
                if len(micro_thresh[0]) == len(df.iloc[0][max_period:min_period]):
                    micro_dead = True

                # Also check if the mean is less than or equal to micro_threshold, indicating a dead channel too,
                # but perhaps not all elements are less than the micro_threshold
                if np.mean(df.iloc[0][max_period:min_period]) <= micro_threshold:
                    micro_dead = True

                ### Trigger on banded character of colorgrid indicating issue ###
                # This could be a dead channel, noise issues, or other issues. This trigger is to let the user know that
                # there's an issue with the health of the station. When we add more complexity to these metrics and
                # review more colorgrids we might be able to get more specific here ###

                # Get rows of dataFrame into a list
                rows = list(df.iloc[0])
                banded_issue = False

                # Then find the average of the difference. Low averages generally indicate that a banded character is
                # present in the colorGrid plot, while high averages are more likely to indicate that the station is
                # healthy
                banded = np.mean(np.diff(rows))

                if banded <= banded_threshold:
                    banded_issue = True

                # Create dictionary for every stream object, then append to data to encompass all information per stream
                d = {
                    "start_time": starttime,
                    "end_time": endtime,
                    "snclq": snclq,
                    "noise_masks": {
                        "frequency_start": freq_start,
                        "frequency_end": freq_end,
                    },
                    "microseism_masks": micro_dead,
                    "banded_masks": banded_issue,
                    "metric_name": "dailyPdfPlot",
                }
                d.update(update_grid)
                d.update(update_line)
                data.append(d)
            except (IndexError, ZeroDivisionError):
                d = {
                    "metric_name": "dailyPdfPlot",
                    "start_time": starttime,
                    "end_time": endtime,
                    "snclq": snclq,
                }
                d.update(update_grid)  # add in the plot image and keys
                d.update(update_line)
                data.append(d)
                print(("dailyPDFPlot Info: No threshold created for SNCLQ: " + snclq))

    if model == "network":
        # Create empty arrays to fill for thresholds
        data = []
        period_masks = []
        freq_start = []
        freq_end = []

        # Create variables for period, enModel, zModel
        if net_dict is not None:
            # Grab out z_power, en_power models and period
            z_power = net_dict["zModel"][1]
            en_power = net_dict["enModel"][1]
            # Get indices of period from the zModel, then convert to a series
            m_period = get_m_period(net_dict)
        else:
            # If a previous dictionary is not defined, create a network noise model on the fly. Recommended that
            # we have many stations in the stream object and > 1 week of data for this

            # Comments inside of networkNoiseModel ask that the network MUST be defined
            # if passing a database instead of stream
            # If this is not included, operations against different Networks will
            # try to be performed/calculations will fail
            net_noise = networkNoiseModel(database, network=network, session=session)
            # Grab out z_power, en_power models and period
            z_power = net_noise["zModel"][1]
            en_power = net_noise["enModel"][1]
            # Get indices of period from the zModel, then convert to a series
            m_period = get_m_period(net_noise)

        tb = database.get_metric(
            "psdMetric",
            network=network,
            station=station,
            channel=channel,
            location=location,
            session=session,
        )

        if tb.empty:
            logger.error(
                "dailyPDFPlot: Error no PSDS found in DB for Network: \
                    {n}, Station: {s}, Channel: {c}, Session: {se}".format(
                    n=network, s=station, c=channel, se=session
                )
            )
            return

        for i in range(len(tb.channel.unique())):
            df_line = pd.DataFrame()
            df = pd.DataFrame()
            date_index = []
            channel = tb.channel.iloc[i]
            tb_psds = tb.uncorrected_psds[tb.channel == channel]
            chan_psds = []
            for j in range(len(tb_psds)):
                psd = database.extract_masks(tb_psds.iloc[j])
                chan_psds.append(psd)
            stats = psdStatistics(chan_psds, logger=logger, database=database)

            for k in range(len(chan_psds)):
                pmode = stats[k]["mode"]
                snclq = tb.snclq[tb.channel == channel].iloc[k]
                date = tb.start_time[tb.channel == channel].iloc[k]
                starttime = tb.start_time[tb.channel == channel].iloc[k]
                endtime = tb.end_time[tb.channel == channel].iloc[k]
                freq = np.asarray(chan_psds[k][0][0])
                period = 1 / freq

                if channel.endswith("Z"):
                    m_power = z_power
                else:
                    m_power = en_power

                date_index.append(date[0:10])
                # loop through all periods and find nearest in gsn model
                diff, per = noiseDiff(period, pmode, m_power, m_period)

                df2 = pd.DataFrame(data=np.asarray(diff).reshape(1, len(diff)), columns=per)
                df = pd.concat([df, df2], ignore_index=True)

                df_line2 = pd.DataFrame(data=pmode.reshape(1, len(pmode)), columns=per)
                df_line = pd.concat([df_line, df_line2], ignore_index=True)

            df["Date"] = date_index
            df = df.rename(df["Date"])
            df = df.drop(["Date"], axis=1)

            df_line["Date"] = date_index
            df_line = df_line.rename(df_line["Date"])
            df_line = df_line.drop(["Date"], axis=1)
            # dropping digitizer values

            df = df.drop([df.columns[-1]], axis=1)
            df = df.drop([df.columns[-2]], axis=1)
            df = df.drop([df.columns[-3]], axis=1)

            update_grid = plot_grid(df, snclq, plot="dailyPDF", f_name=f_name_grid)
            update_line = _plot_colorline(df_line, snclq, per_arr, f_name=f_name_line)

            try:
                ### Trigger based on large differences from noise model ###
                # Trigger mask output based on diff_threshold value -- looks for
                thresh = np.where(df >= diff_threshold)[1]
                # Get consecutive indices out
                thresh_cons = consecutive(thresh, stepsize=1)

                # Loop through consecutive indices and find where freq values where threshold exceeds those values
                for k in range(len(thresh_cons)):
                    period_masks.append(df.columns[thresh_cons[k]])

                # Grab out period masks and aggregate them into start/stop frequencies
                for k in range(len(period_masks)):
                    freq_start.append(1 / float(period_masks[k][0]))
                    freq_end.append(1 / float(period_masks[k][-1]))

                ### Trigger based on unhealthy differences from noise model in the microseism band ###
                # Microseism band is 0.1 - 1.0 Hz, which translates to 1 - 10s periods
                # Find indices that are between 1-10s period so that we can subset the dataframe accordingly
                # First convert columns to list, then convert strings to floats
                cols = list(df.columns)
                cols_f = [float(i) for i in cols]
                micro_dead = False
                # Find index of max period (10s) and min period (1s)
                max_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 10))
                min_period = min(list(range(len(cols_f))), key=lambda i: abs(cols_f[i] - 1))

                # Find where values are below threshold, indicating a potentially dead channel
                micro_thresh = np.where(df.iloc[0][max_period:min_period] <= micro_threshold)

                # Check if length of micro_thresh equal to number of elements in period range of microseism, this
                # means that all of them are less than or equal to the threshold, indicating a dead channel
                if len(micro_thresh[0]) == len(df.iloc[0][max_period:min_period]):
                    micro_dead = True

                # Also check if the mean is less than or equal to micro_threshold, indicating a dead channel too,
                # but perhaps not all elements are less than the micro_threshold
                if np.mean(df.iloc[0][max_period:min_period]) <= micro_threshold:
                    micro_dead = True

                ### Trigger on banded character of colorgrid indicating issue ###
                # This could be a dead channel, noise issues, or other issues. This trigger is to let the user know that
                # there's an issue with the health of the station. When we add more complexity to these metrics and
                # review more colorgrids we might be able to get more specific here ###

                # Get rows of dataFrame into a list
                rows = list(df.iloc[0])
                banded_issue = False

                # Then find the average of the difference. Low averages generally indicate that a banded character
                # is present in the colorGrid plot, while high averages are more likely to indicate that the station
                # is healthy
                banded = np.mean(np.diff(rows))

                if banded <= banded_threshold:
                    banded_issue = True

                # Create dictionary for every stream object, then append to data to encompass all information per
                # stream
                d = {
                    "start_time": starttime,
                    "end_time": endtime,
                    "snclq": snclq,
                    "noise_masks": {
                        "frequency_start": freq_start,
                        "frequency_end": freq_end,
                    },
                    "microseism_masks": micro_dead,
                    "banded_masks": banded_issue,
                    "metric_name": "dailyPdfPlot",
                }
                d.update(update_grid)
                d.update(update_line)
                data.append(d)

            except (IndexError, ZeroDivisionError):
                d = {
                    "metric_name": "dailyPdfPlot",
                    "start_time": starttime,
                    "end_time": endtime,
                    "snclq": snclq,
                }
                d.update(update_grid)  # add in the plot image and keys
                d.update(update_line)
                data.append(d)
                print(("dailyPDFPlot Info: No threshold created for SNCLQ: " + snclq))
                continue
    database.insert_metric(data)
    return data


def _getDays(tr):
    """
    Splits a multi-day trace into individual days. Each day is returned as a new trace.

    :param tr: trace objects
    :type tr: obspy.core.trace.Trace

    :return: An obspy Stream object containing single day traces.
    :rtype: obspy.core.stream.Stream

    """
    x = tr.data
    samp = tr.stats.sampling_rate
    day = samp * (60 * 60 * 24)  # number of samples in 1 day
    sep_data = []
    count = 0
    # looping through samples by number of samples in a day
    while count < len(x):
        # spliting data into 1day samples
        split = x[int(count) : int(count + day)]
        sep_data.append(split)
        # adding samples in day to count to increment
        count = int(count + day)
    # creating new stream
    stN = obspy.Stream()
    # looping through separated data
    for i in range(len(sep_data)):
        # making sure not to change original starttime of first trace
        if i != 0:
            tr.stats["starttime"] = tr.stats.starttime + (60 * 60 * 24)
        tr.stats["npts"] = len(sep_data[i])
        # creating new trace
        trN = obspy.Trace(sep_data[i], header=tr.stats)
        # adding new trace to stream
        stN.append(trN)
    # making sure that each traces has the same number of samples (each 1 full day),
    #  otherwise psdMetric will throw error
    for j in range(1, len(stN)):
        trN1 = stN[j]
        trN2 = stN[j - 1]
        if trN1.stats["npts"] != trN2.stats["npts"]:
            stN.remove(trN1)
    return stN


def _find_nearest(array, value):
    """
    Finds the closest value in a given array

    :param array: Numpy Array containing  numbers
    :type array: numpy.ndarray
    :param value: The value to find the closest to.
    :type value: float or int

    :return: The index of the closest value in the array
    :rtype: int


    """
    if array.dtype != float:
        array = np.array(array, dtype=float)
    idx = (np.abs(array - float(value))).argmin()
    return idx


def noiseDiff(period, mode, m_power, m_period):
    """
    Calculates differences from model power mode and computed power mode

    :param period:  Computed period
    :type period: numpy.ndarray
    :param mode: Computed mode
    :type mode: numpy.ndarray
    :param m_power: Model power
    :type m_power: numpy.ndarray
    :param m_period: Model period
    :type m_period: numpy.ndarray
    :return: a list of periods and a list of differences

    """
    # loop through all periods and find nearest in gsn model
    try:
        diff = []
        per = []
        for k in range(len(period)):
            index = _find_nearest(np.asarray(m_period), period[k])
            diff.append(mode[k] - m_power[index])
            per_str = str(period[k])
            per.append(per_str[0:5])
        return diff, per
    except IndexError:
        diff = []
        per = []
        for k in range(len(mode)):
            index = _find_nearest(np.asarray(m_period), period[k])
            diff.append(mode[k] - m_power[index])
            per_str = str(period[k])
            per.append(per_str[0:5])
        return diff, per


def plot_grid_data_fill_in(grd_data):
    # Check to make sure all rows have Unique names (possibly by channel?)
    plt.style.use("ggplot")
    # Binning data into color ranges
    grd_data = grd_data.fillna(99999)

    for i in grd_data.columns:
        for j in grd_data.index:
            try:
                if grd_data.at[j, i] == 99999.0:
                    grd_data.at[j, i] = 80
            except ValueError:
                grd_data = grd_data.drop_duplicates()
            if grd_data.at[j, i] == 99999:
                grd_data.at[j, i] = 80

            elif grd_data.at[j, i] >= 70 and grd_data.at[j, i] < 99999:
                grd_data.at[j, i] = 70

            elif grd_data.at[j, i] >= 60 and grd_data.at[j, i] < 70:
                grd_data.at[j, i] = 60

            elif grd_data.at[j, i] >= 50 and grd_data.at[j, i] < 60:
                grd_data.at[j, i] = 50

            elif grd_data.at[j, i] >= 40 and grd_data.at[j, i] < 50:
                grd_data.at[j, i] = 40

            elif grd_data.at[j, i] >= 30 and grd_data.at[j, i] < 40:
                grd_data.at[j, i] = 30

            elif grd_data.at[j, i] >= 20 and grd_data.at[j, i] < 30:
                grd_data.at[j, i] = 20

            elif grd_data.at[j, i] >= 10 and grd_data.at[j, i] < 20:
                grd_data.at[j, i] = 10

            elif grd_data.at[j, i] >= 0 and grd_data.at[j, i] < 10:
                grd_data.at[j, i] = 0

            else:
                grd_data.at[j, i] = -10

    # flipping data so that period goes from small to large, and date newer to older
    grd_data = grd_data.reindex(columns=grd_data.columns[::-1])
    return grd_data


def plot_grid(grd_data, snclq, plot, f_name=None, channel=None, rank=None):
    """
    Plots colorgrid data into Colorgrid (power mode or station ranking). Colorgrid plots show the magnitude of deviation
    (in decibels (dB)) of a station's daily PSD PDF mode from the New Low Noise Model (NLNM) (Peterson, 1993), GSN 2004
    model, or a user provided/calculated network noise model. Station ranking plots show the magnitude of deviation
    (in dB) of the mean of each station’s daily PDF mode  (averaged at each period bin over the respective time period)
    from the NLNM. Stations are ranked from quietest to noisiest based on user provided frequencies.  Hotter colors
    indicate noisier data, while cooler colors indicate quieter data on colorgrid and station ranking plots.

    :param grd_data: (dataframe) dataframe containing period and values
    :param snclq: (str) SNCLQ id
    :param plot: (str) Plot type (stationRanking or dailyPDF)
    :param f_name: (str) plot filename to save to. If None, image will not be saved.
    :param channel: (str) Channel name
    :param rank: (arr) Array containing ranking values.

    :return: Station Ranking or Daily PDF Colorgrid plot

    Colorgrid plots, daily PDF mode timeline plots, and station ranking plots were first established for use in the
    USArray TA project (Busby et. al., 2018) and code was ported from Perl and Hypertext Preprocessor (PHP) to Python.
    Code was adapted and augmented for use within Pycheron.

    #. Busby, R. W., R.L. Woodward, K.A. Hafner, F. L. Vernon, and A.M. Frasseto (2018). The Design and Implementation
    of EarthScope's USArray Transportable Array in the Conterminous United States and Southern Canada (Rep.).

    """
    plt.style.use("ggplot")
    # Binning data into color ranges
    grd_data = grd_data.fillna(99999)

    for i in grd_data.columns:
        for j in grd_data.index:
            try:
                if grd_data.at[j, i] == 99999.0:
                    grd_data.at[j, i] = 80
            except ValueError:
                grd_data = grd_data.drop_duplicates()

            if grd_data.at[j, i] == 99999:
                grd_data.at[j, i] = 80

            elif grd_data.at[j, i] >= 70 and grd_data.at[j, i] < 99999:
                grd_data.at[j, i] = 70

            elif grd_data.at[j, i] >= 60 and grd_data.at[j, i] < 70:
                grd_data.at[j, i] = 60

            elif grd_data.at[j, i] >= 50 and grd_data.at[j, i] < 60:
                grd_data.at[j, i] = 50

            elif grd_data.at[j, i] >= 40 and grd_data.at[j, i] < 50:
                grd_data.at[j, i] = 40

            elif grd_data.at[j, i] >= 30 and grd_data.at[j, i] < 40:
                grd_data.at[j, i] = 30

            elif grd_data.at[j, i] >= 20 and grd_data.at[j, i] < 30:
                grd_data.at[j, i] = 20

            elif grd_data.at[j, i] >= 10 and grd_data.at[j, i] < 20:
                grd_data.at[j, i] = 10

            elif grd_data.at[j, i] >= 0 and grd_data.at[j, i] < 10:
                grd_data.at[j, i] = 0

            else:
                grd_data.at[j, i] = -10

    # flipping data so that period goes from small to large, and date newer to older
    grd_data = grd_data.reindex(columns=grd_data.columns[::-1])

    if plot == "stationRanking":
        grd_data = grd_data.reindex(index=grd_data.index[::-1])
        # initializing plot
        if len(grd_data) < 10:
            fig, ax = plt.subplots(figsize=(10, 5))
        else:
            fig, ax = plt.subplots(figsize=(10, (len(grd_data) * 0.25)))
    else:
        if len(grd_data) < 15:
            # initializing plot
            fig, ax = plt.subplots(figsize=(10, 5))
        else:
            # initializing plot
            fig, ax = plt.subplots(figsize=(10, (len(grd_data) * 0.2)))

    # making custom colorbar that is discrete
    cmap = matplotlib.colors.ListedColormap(
        [
            "#0000FF",
            "#3366FF",
            "#66CCFF",
            "#66FFCC",
            "#CCFF66",
            "orange",
            "red",
            "#993366",
            "#660033",
            "black",
        ]
    )
    # cmaplist = [cmap(i) for i in range(cmap.N)]
    # cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = np.linspace(0, 10, 11)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    # Plotting Data
    ax.pcolor(grd_data, cmap=cmap, vmin=-10, vmax=80)
    # Format
    fig = plt.gcf()

    ax.set_aspect("equal", adjustable="box")

    # Crfeating colorbar
    ax2 = fig.add_axes([0.110, 0.84, 0.8, 0.025])
    matplotlib.colorbar.ColorbarBase(
        ax2,
        cmap=cmap,
        norm=norm,
        spacing="proportional",
        ticks=np.linspace(0.5, 9.5, 10),
        boundaries=bounds,
        format="%1i",
        orientation="horizontal",
    )
    ax2.set_xticklabels(
        [
            "D<0",
            "10>D>=0",
            "20>D>=10",
            "30>D>=20",
            "40>D>=30",
            "50>D>=40",
            "60>D>=50",
            "70>D>=60",
            "D>=70",
            "No Data",
        ],
        size=7,
    )

    ax.set_xlabel("Period (seconds)", size=12)
    try:
        snclq_id = snclq.split(".")
    except AttributeError:
        snclq_id = snclq[0].split(".")

    if plot == "dailyPDF":
        dpi = 300
        ax.set_ylabel("Date", size=12)
        plt.title(
            "Daily PDF Mode Power Grid: " + snclq_id[0] + " " + snclq_id[1] + " " + snclq_id[2] + " " + snclq_id[3],
            y=4,
            weight="bold",
            size=12,
        )

        # adding axis labels and ticks
        ax.yaxis.set_ticks(np.arange(0.5, len(grd_data.index) + 0.5, 1))
        ax.set_yticklabels(grd_data.index, size=4)

    if plot == "stationRanking":
        ax.set_ylabel("Station", size=12)

        plt.title(
            "Station Ranking Plot for Network: " + snclq_id[0] + ", Ranked by: " + str(rank) + "s",
            y=4,
            weight="bold",
            size=12,
        )
        fig.suptitle(
            "                      Channel: " + channel,
            y=0.94,
            size=10,
            horizontalalignment="center",
        )
        ax.yaxis.set_ticks(np.arange(0.5, len(grd_data.index) + 0.5))
        if len(grd_data.index) > 100:
            dpi = 400
            ax.set_yticklabels(grd_data.index, size=4)
        else:
            dpi = 300
            ax.set_yticklabels(grd_data.index, size=7)

        # draw box around ranked cell
        x = list(grd_data.columns).index(str(rank))
        y = len(grd_data.index)
        w = 1

        points = [[x, 0], [x + w, 0], [x + w, y], [x, y]]
        box = plt.Polygon(points, fill=None, snap=True, lw=1.5, ec="k")
        ax.add_patch(box)

    ax.xaxis.set_ticks(np.arange(0.5, len(grd_data.columns), 2))
    ax.set_xticklabels(grd_data.columns[0::2], size=7, rotation="vertical")

    ax.tick_params(
        axis="both",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,
        left=False,  # ticks along the bottom edge are off
        right=False,
        labeltop=True,
    )  # ticks along the top edge are off
    ax2.tick_params(
        axis="both",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,
        left=False,  # ticks along the top edge are off
        labelbottom=False,
        labeltop=True,
        labelleft=False,
    )

    if f_name is not None and plot == "dailyPDF":
        name = f_name + "_" + snclq_id[3] + ".png"
        plt.savefig(name, dpi=dpi)
        d = {"pdf_grid_image_path": name}
        return d

    elif f_name is not None and plot == "stationRanking":
        name = f_name + "_" + snclq_id[0] + ".png"
        plt.savefig(name, dpi=dpi)
        d = {"pdf_grid_image_path": name}
        return d
    else:
        plt.show()


def _plot_colorline(grd_data, snclq, per_arr=None, f_name=None):
    """

    Creates a daily PDF mode timeline plot. PDF mode timeline plots show the daily mode power at a set of user defined
    frequencies over a specified time interval. These plots are utilized to track background noise relative to the
    frequencies of interest.

    :param grd_data: (dataframe) Pandas dataframe of data.
    :param snclq:(str) sncql id
    :param per (array): Array of periods to plot
    :param f_name (str): file name plot png

    Colorgrid plots, daily PDF mode timeline plots, and station ranking plots were first established for use in the
    USArray TA project (Busby et. al., 2018) and code was ported from Perl and Hypertext Preprocessor (PHP) to Python.
    Code was adapted and augmented for use within Pycheron.

    #. Busby, R. W., R.L. Woodward, K.A. Hafner, F. L. Vernon, and A.M. Frasseto (2018). The Design and Implementation
    of EarthScope's USArray Transportable Array in the Conterminous United States and Southern Canada (Rep.).

    """

    if per_arr is None:
        # default values
        per = [0.1, 1.0, 6.5, 10, 30.8, 103.7]
    else:
        per = per_arr

    plt.style.use("ggplot")
    fig, ax = plt.subplots(figsize=(11, 8))
    # Format
    fig = plt.gcf()
    fig.set_size_inches(11, 8)
    for i in per:
        freq_index = _find_nearest(np.asarray(grd_data.columns), i)
        ax.plot(
            grd_data[grd_data.columns[freq_index]],
            label=str(grd_data.columns[freq_index])[0:5] + " Sec",
        )
    plt.legend(
        bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
        loc=3,
        ncol=len(per),
        mode="expand",
        borderaxespad=0.0,
    )

    ax.xaxis.set_ticks(np.arange(0.5, len(grd_data.index), 2))
    ax.set_xticklabels(grd_data.index[0::2], size=7, rotation="vertical")

    plt.subplots_adjust(bottom=0.2, top=0.85)
    plt.ylabel("Power (dB)")
    plt.xlabel("Date")
    snclq_id = snclq.split(".")
    plt.title(
        "Daily PDF Mode Power Timelines: " + snclq_id[0] + " " + snclq_id[1] + " " + snclq_id[2] + " " + snclq_id[3],
        y=1.15,
    )

    if f_name is not None:
        name = f_name + "_" + snclq_id[3] + ".png"
        plt.savefig(name)
        d = {"pdf_line_image_path": name}
        return d
