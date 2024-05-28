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

from pycheron.db.sqllite_db import Database
from pycheron.psd.psdList import psdList
from pycheron.psd.psdStatistics import psdStatistics
from pycheron.util.logger import Logger
from pycheron.util.format import parse_snclq
from pycheron.util.getChannel import getChannelName


plt.style.use("ggplot")

__all__ = ["stationNoiseModel"]


def stationNoiseModel(
    st,
    type="envelope10_90",
    plot=False,
    fname=None,
    network=None,
    station=None,
    location=None,
    session=None,
    database=None,
    logger=None,
):
    """
    Calculates Station Noise Model for a station

    :param st: Stream object of data from one station or a Database object
    :type st: obspy.core.stream.Stream or pycheron.db.sqllite_db.Database
    :param type: Type of model:

        * "single10" : Based on single line of 10th percentile
        * "single90" : Based on single line of 90th percentile
        * "single05" : Based on single line of 5th percentile
        * "single95" : Based on single line of 95th percentile
        * "envelope10_90" : Based on envelope between 10th and 90th percentile
        * "envelope05_95" : Based on envelope between 5th and 95th percentile

    :type type: str
    :param plot: Whether to output plot
    :type plot: bool
    :param fname: If plot = True, the file location to save it
    :type fname: str
    :param network: If using a db object, the network of the data you want to use
    :type network: str
    :param station: If using a db object, the station of the data you want to use
    :type station: str
    :param location: If using a db object, the location of the data you want to use
    :type location: str
    :param session: If using a db object, the session name of the data you want to use
    :type session: str
    :param database:  database object to save to
    :type database: pycheron.db.sqllite_db.Database
    :param logger: Logger object
    :type logger: pycheron.util.logger.Logger

    :return: A dictionary with the following keys and types:

        * metric_name (`str`)
        * network (`str`)
        * station (`str`)
        * location (`str`)
        * type (`str`)
        * If *type* is: "single10", "single90", "single05", or "single95":

          * "e_station_noiseModel" (`pandas.Series`) - station noise model based on percentile for channel \*\*E
          * "n_station_noiseModel" (`pandas.Series`) - station noise model based on percentile for channel \*\*N
          * "z_station_noiseModel" (`pandas.Series`) - station noise model based on percentile for channel \*\*Z

        * If *type* is: "envelope10_90" or "envelope05_95":

          * "e_station_noiseModel_low" (`pandas.Series`) - station noise model based on higher bound percentile
            (90 or 95) for channel \*\*E
          * "e_station_noiseModel_high" (`pandas.Series`) -station noise model based on lower bound percentile
            (05 or 10) for channel \*\*E
          * "n_station_noiseModel_low" (`pandas.Series`) - station noise model based on higher bound percentile
            (90 or 95) for channel \*\*N
          * "n_station_noiseModel_high" (`pandas.Series`) -station noise model based on lower bound percentile
            (05 or 10) for channel \*\*N
          * "z_station_noiseModel_low" (`pandas.Series`) - station noise model based on higher bound percentile
            (90 or 95) for channel \*\*Z
          * "z_station_noiseModel_high" (`pandas.Series`) -station noise model based on lower bound percentile
            (05 or 10) for channel \*\*Z

    :rtype: dict

    Station Noise models are based on a single line (based off a statistic) or an envelope. The envelope method follows
    the method outlined in this paper: [#]_

    ** Algorithm Steps **

    #. Compute 5th, 10th, 90th, 95th, and 50th percentiles (part of psdStatistics function)
    #. The station noise baseline "envelope" is basd off the 10th and 90th percentiles of the PSD PDF (Provides options
       to include 5th, 10th, 90th, 95th, and an option for envelope vs. single
    #. For the single line option, do it based on the mode at all periods of the station, following the idea of the
       network noise model, so it's the mean of the modes. All channel components have their own model.

    .. rubric:: References
       .. [#] McNamara, D. E., Hutt, C. R., Gee, L. S., Benz, H. M., & Buland, R. P. (2009). A method to establish
              seismic noise baselines for automated station assessment. Seismological Research Letters, 80(4), 628-637.

    **Example**

    .. code-block:: python

        from pycheron.psd.noise.stationNoiseModel import stationNoiseModel
        from obspy import UTCDateTime
        from obspy.clients.fdsn import Client


        client = Client("IRIS")
        t = UTCDateTime("2017-06-01T00:00:00")
        te = UTCDateTime("2017-06-02T00:00:00")
        st = client.get_waveforms("TX", "ALPN", "00", "*", t, te)

        model = stationNoiseModel(st, type="envelope10_90",plot=True)

        # access z station noise above 90 percent
        model['z_station_noise_high']
        # access z station noise below 10 percent
        model['z_station_noise_low']
        # access e station noise above 90 percent
        model['e_station_noise_high']
        # access e station noise below 10 percent
        model['e_station_noise_low']
        # access n station noise above 90 percent
        model['n_station_noise_high']
        # access n station noise below 10 percent
        model['n_station_noise_low']

    **Plotting**

    If you set the parameter `plot` equal to `True` then the following plot will be produced.

    .. code-block:: python

        # using output `model` from above example
        import matplotlib.pyplot as plt

        # Station Noise Plots
        plt.semilogx(model["e_station_noise_low"], 'orangered', label="E Station")
        plt.semilogx(model["e_station_noise_high"], 'orangered')

        plt.semilogx(model["n_station_noise_low"], 'c', label="N Station")
        plt.semilogx(model["n_station_noise_high"], 'c')

        plt.semilogx(model["z_station_noise_low"], 'darkmagenta', label="Z Station")
        plt.semilogx(model["z_station_noise_high"], 'darkmagenta')

        plt.xlabel("Period")
        plt.ylabel("Power(dB)")
        plt.title("Noise Model for Station: " + st[0].stats.station + " (" + "envelope10_90" + ")")

        plt.grid(True, "both", "both", c="white")
        plt.legend()

    .. image:: _static/stationNoiseModel.png
    """
    # import pdb; pdb.set_trace();

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # If obspy stream object, calculate psds and grab out snclq info
    if isinstance(st, obspy.core.stream.Stream):
        network, station, channel, location, quality = parse_snclq(st[0].get_id())

        # If database object, check to see if we already have the psds calculated
        if database is not None:
            tb = database.get_metric(
                "psdMetric",
                network=network,
                station=station,
                channel=channel,
                session=session,
            )

            # if returns no data, calculate and insert into db
            if tb.empty:
                print(f"TB IS EMPTY")
                psds = psdList(st)
                print(f"@@@@@@@@ PSDS: {psds}")
                database.insert_metric(psds)

            else:
                print(f"TB NOT EMPTY")
                tb_psds = tb.uncorrected_psds[tb.channel == channel]
                print(f"TB PSDS: {tb_psds}")
                psds = []
                for j in range(len(tb_psds)):
                    psd = database.extract_masks(tb_psds.iloc[j])
                    psds.append(psd)
                print(f"########## PSDS: {psds}")
        # Otherwise, calculate psds
        else:
            psds = psdList(st)
            print(f"$$$$$$$$$$$$ PSDS: {psds}")

    # If st is a database object, then grab psd from db. Grab out station name.
    else:
        database = st
        psds = create_psds_from_metric_table(database, network, station, session, location, logger)
        print(f"****************** PSDS: {psds}")

    # get psd statistics for the psds
    print(f"!!!!!!!! PSDS: {psds}")
    df_dict = gen_stats_plot_vals(psds, type, logger, database)

    # Based on theose statistics, generate single or envelope list
    out = simple_or_envelope(df_dict, network, station, type, location)

    # Add dfs created to output dictionary. If type startswith "s", meaning we are doing single noise model lines,
    # add to output dictionary and plot if plot == True, plot the single profiles using orange red for the E component,
    # cyan for the N component and dark magenta for the Z component. Make axis labels, grid, and legend and save figure
    # to fname
    if type.startswith("s"):
        if plot:
            f = plt.Figure()
            p = f.add_subplot(1, 1, 1)

            # Station Noise Plots
            p.semilogx(out["e_station_noiseModel"], "orangered", label="E Station")

            p.semilogx(out["n_station_noiseModel"], "c", label="N Station")

            p.semilogx(out["z_station_noiseModel"], "darkmagenta", label="Z Station")

            p.set_xlabel("Period")
            p.set_ylabel("Power(dB)")
            p.set_title("Noise Model for Station: " + station + " (" + type + ")")
            p.set_xlim(0.01, 200)
            p.grid(True, "both", "both", c="white")
            p.legend()
            if fname is not None:
                f.savefig(fname)
                update = {"station_noise_image_path": fname}
                out.update(update)
                plt.close()
        if isinstance(st, Database) or database is not None:
            database.insert_metric(out)
        else:
            return out

    # Add dfs created to output dictionary. If type doesn't start with "s", means we are doing envelopes so plot
    # high/low noise profiles, add to output dictionary and plot if plot == True, plot the envelopes using orange red
    # for the E component high/low noise profile, cyan for the N component high/low noise profile, and dark magenta for
    # the Z component high/low noise profile. Make axis labels, grid, and legend and save figure to
    # fname
    else:
        if plot:
            f = plt.figure()
            p = f.add_subplot(1, 1, 1)

            # Station Noise Plots
            p.semilogx(out["e_station_noiseModel_low"], "orangered", label="E Station")
            p.semilogx(out["e_station_noiseModel_high"], "orangered")

            p.semilogx(out["n_station_noiseModel_low"], "c", label="N Station")
            p.semilogx(out["n_station_noiseModel_high"], "c")

            p.semilogx(out["z_station_noiseModel_low"], "darkmagenta", label="Z Station")
            p.semilogx(out["z_station_noiseModel_high"], "darkmagenta")

            p.set_xlabel("Period")
            p.set_ylabel("Power(dB)")
            p.set_title("Noise Model for Station: " + station + " (" + type + ")")
            p.set_xlim(0.01, 200)
            p.grid(True, "both", "both", c="white")
            p.legend()

            if fname is not None:
                f.savefig(fname)
                update = {"station_noise_image_path": fname}
                out.update(update)
                plt.close()
        if isinstance(st, Database) or database is not None:
            database.insert_metric(out)
        else:
            return out


def create_psds_from_metric_table(database, network, station, session=None, location=None, logger=None):
    """
    Grab psds from database
    """
    # Grab psdMetric information out of the databsae
    tb = database.get_metric(
        "psdMetric",
        network=network,
        station=station,
        session=session,
        channel=None,
        location=location,
    )
    print(f"TB: {tb}")
    # If empty, complain and exit
    if tb.empty:
        logger.error(
            "stationNoiseModel Error: No psdMetric found for specific network in Database, please run psdMetric"
            " first"
        )
        return
    # If location null, add dashes
    if location is None:
        location = "--"

    # Initialize empty psd list
    psds = []
    print(f"Length of Uncorrected_psds: {tb.uncorrected_psds}")
    # Loop through psds an append to new output list
    for i in range(len(tb.uncorrected_psds)):
        psd = database.extract_masks(tb.uncorrected_psds.iloc[i])
        psds.append(psd)

    return psds


def gen_stats_plot_vals(psds, type, logger=None, database=None):
    """
    Generate psd statistics values for station network noise plots
    """
    # initialize dataframes to fill
    df_z_en1 = pd.DataFrame()
    df_z_en2 = pd.DataFrame()
    df_e_en1 = pd.DataFrame()
    df_e_en2 = pd.DataFrame()
    df_n_en1 = pd.DataFrame()
    df_n_en2 = pd.DataFrame()
    # initializing single dfs to fill
    df_e = pd.DataFrame()
    df_n = pd.DataFrame()
    df_z = pd.DataFrame()

    # Calculate psd statistics
    stats = psdStatistics(psds, logger=logger, database=database)

    # Try grabbing out period
    try:
        per = 1 / psds[0][0][0]
    except IndexError:
        return
    except TypeError:
        per = 1 / np.asarray(psds[0][0][0])

    # looping through each psd
    for i in range(len(stats)):

        chan = getChannelName(stats[i]["snclq"])
        # makes sure no nan values in stats
        if (
            not np.isnan(stats[i]["percent_5"]).all()
            and not np.isnan(stats[i]["percent_10"]).all()
            and not np.isnan(stats[i]["percent_90"]).all()
            and not np.isnan(stats[i]["percent_95"]).all()
        ):
            # Get percentages out
            per05 = stats[i]["percent_5"]
            per10 = stats[i]["percent_10"]
            per90 = stats[i]["percent_90"]
            per95 = stats[i]["percent_95"]

            ########
            # Do we need this code too? If we know that period is the same length as stats[i]["*Percent"] do we need
            # this? Or are there situations we found where they aren't?
            #######
            # if len(per) <= len(stats[i]["percent_5"]):
            #     per05 = stats[i]["percent_5"][np.where(per)]
            #     per10 = stats[i]["percent_10"][np.where(per)]
            #     per90 = stats[i]["percent_90"][np.where(per)]
            #     per95 = stats[i]["percent_95"][np.where(per)]
            # else:
            #     per = per[np.where(stats[i]["percent_5"])]
            #     per05 = per
            #     per10 = per
            #     per90 = per
            #     per95 = per

            # For channels that endwith Z, structure the dataframes appropriately for the various options: single10,
            # single05, single90, single95, or the envelopes of 10 and 90 percentiles or 5 and 95 percentiles
            if chan.endswith("Z"):
                if type == "single10":
                    df2 = pd.DataFrame(data=per10.reshape(1, len(per10)), columns=per)
                    df_z = pd.concat([df_z, df2], ignore_index=True)
                elif type == "single05":
                    df2 = pd.DataFrame(data=per05.reshape(1, len(per05)), columns=per)
                    df_z = pd.concat([df_z, df2], ignore_index=True)
                elif type == "single90":
                    df2 = pd.DataFrame(data=per90.reshape(1, len(per90)), columns=per)
                    df_z = pd.concat([df_z, df2], ignore_index=True)
                elif type == "single95":
                    df2 = pd.DataFrame(data=per95.reshape(1, len(per95)), columns=per)
                    df_z = pd.concat([df_z, df2], ignore_index=True)
                elif type == "envelope10_90":
                    df2 = pd.DataFrame(data=per10.reshape(1, len(per10)), columns=per)
                    df_z_en1 = pd.concat([df_z_en1, df2], ignore_index=True)
                    df2 = pd.DataFrame(data=per90.reshape(1, len(per90)), columns=per)
                    df_z_en2 = pd.concat([df_z_en2, df2], ignore_index=True)
                elif type == "envelope05_95":
                    df2 = pd.DataFrame(data=per05.reshape(1, len(per05)), columns=per)
                    df_z_en1 = pd.concat([df_z_en1, df2], ignore_index=True)
                    df2 = pd.DataFrame(data=per95.reshape(1, len(per95)), columns=per)
                    df_z_en2 = pd.concat([df_z_en2, df2], ignore_index=True)

            # For channels that endwith E or 2, structure the dataframes appropriately for the various options:
            # single10, single05, single90, single95, or the envelopes of 10 and 90 percentiles or 5 and 95 percentiles
            elif chan.endswith("E") or chan.endswith("2"):

                if type == "single10":
                    df2 = pd.DataFrame(data=per10.reshape(1, len(per10)), columns=per)
                    df_e = pd.concat([df_e, df2], ignore_index=True)
                elif type == "single05":
                    df2 = pd.DataFrame(data=per05.reshape(1, len(per05)), columns=per)
                    df_e = pd.concat([df_e, df2], ignore_index=True)
                elif type == "single90":
                    df2 = pd.DataFrame(data=per90.reshape(1, len(per90)), columns=per)
                    df_e = pd.concat([df_e, df2], ignore_index=True)
                elif type == "single95":
                    df2 = pd.DataFrame(data=per95.reshape(1, len(per95)), columns=per)
                    df_e = pd.concat([df_e, df2], ignore_index=True)
                elif type == "envelope10_90":
                    df2 = pd.DataFrame(data=per10.reshape(1, len(per10)), columns=per)
                    df_e_en1 = pd.concat([df_e_en1, df2], ignore_index=True)
                    df2 = pd.DataFrame(data=per90.reshape(1, len(per90)), columns=per)
                    df_e_en2 = pd.concat([df_e_en2 ,df2], ignore_index=True)
                elif type == "envelope05_95":
                    df2 = pd.DataFrame(data=per05.reshape(1, len(per05)), columns=per)
                    df_e_en1 = pd.concat([df_e_en1, df2], ignore_index=True)
                    df2 = pd.DataFrame(data=per95.reshape(1, len(per95)), columns=per)
                    df_e_en2 = pd.concat([df_e_en2 ,df2], ignore_index=True)

            # For channels that endwith N or 1, structure the dataframes appropriately for the various options:
            # single10, single05, single90, single95, or the envelopes of 10 and 90 percentiles or 5 and 95 percentiles
            elif chan.endswith("N") or chan.endswith("1"):
                if type == "single10":
                    df2 = pd.DataFrame(data=per10.reshape(1, len(per10)), columns=per)
                    df_n = pd.concat([df_n, df2], ignore_index=True)
                elif type == "single05":
                    df2 = pd.DataFrame(data=per05.reshape(1, len(per05)), columns=per)
                    df_n = pd.concat([df_n, df2], ignore_index=True)
                elif type == "single90":
                    df2 = pd.DataFrame(data=per90.reshape(1, len(per90)), columns=per)
                    df_n = pd.concat([df_n, df2], ignore_index=True)
                elif type == "single95":
                    df2 = pd.DataFrame(data=per95.reshape(1, len(per95)), columns=per)
                    df_n = pd.concat([df_n, df2], ignore_index=True)
                elif type == "envelope10_90":
                    df2 = pd.DataFrame(data=per10.reshape(1, len(per10)), columns=per)
                    df_n_en1 = pd.concat([df_n_en1, df2], ignore_index=True)
                    df2 = pd.DataFrame(data=per90.reshape(1, len(per90)), columns=per)
                    df_n_en2 = pd.concat([df_n_en2, df2], ignore_index=True)
                elif type == "envelope05_95":
                    df2 = pd.DataFrame(data=per05.reshape(1, len(per05)), columns=per)
                    df_n_en1 = pd.concat([df_n_en1, df2], ignore_index=True)
                    df2 = pd.DataFrame(data=per95.reshape(1, len(per95)), columns=per)
                    df_n_en2 = pd.concat([df_n_en2, df2], ignore_index=True)

    # Create output dictionary of values
    df_dict = {}
    df_dict["df_z"] = df_z
    df_dict["df_z_en1"] = df_z_en1
    df_dict["df_z_en2"] = df_z_en2
    df_dict["df_e"] = df_e
    df_dict["df_e_en1"] = df_e_en1
    df_dict["df_e_en2"] = df_e_en2
    df_dict["df_n"] = df_n
    df_dict["df_n_en1"] = df_n_en1
    df_dict["df_n_en2"] = df_n_en2

    return df_dict


def simple_or_envelope(df_dict, network, station, type, location=None):
    """
    Depending out type specified output dictionary for single network noise model line or envelope
    """

    # Single output dictionary
    if type.startswith("s"):
        out = {
            "e_station_noiseModel": df_dict["df_e"].mean(),
            "n_station_noiseModel": df_dict["df_n"].mean(),
            "z_station_noiseModel": df_dict["df_z"].mean(),
            "metric_name": "stationNoiseModel",
            "network": network,
            "station": station,
            "location": location,
            "type": type,
        }
    # Envelope output dictionary
    else:
        out = {
            "e_station_noiseModel_low": df_dict["df_e_en1"].mean(),
            "e_station_noiseModel_high": df_dict["df_e_en2"].mean(),
            "n_station_noiseModel_low": df_dict["df_n_en1"].mean(),
            "n_station_noiseModel_high": df_dict["df_n_en2"].mean(),
            "z_station_noiseModel_low": df_dict["df_z_en1"].mean(),
            "z_station_noiseModel_high": df_dict["df_z_en2"].mean(),
            "metric_name": "stationNoiseModel",
            "network": network,
            "station": station,
            "location": location,
            "type": type,
        }
    return out
