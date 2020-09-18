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

"""
Created on Tue Oct 4 14:58:00 2016
@author: kaaur

PSD plot, useful for plotting instrument corrrected noise values of PSDs and 
PDFs from pycheron package. Creates visualizations for sets of PSDs. Plots generated 
with style ='pdf' mimic plots presented in McNamara paper. Requires data have 
already been read in

Modified by jbobeck (3.5.18)
    -Commented out plt.hold(True), as it is deprecated and giving warnings.
"""

__all__ = ["psdPlot", "_psdPlot_from_database", "_read_psds", "get_psd_plot_data"]

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pycheron.psd.psdStatistics import psdStatistics
from pycheron.psd.psdList import psdList
from pycheron.psd.noise.getNoise import getNoise
from pycheron.psd.noise.noiseModel import noiseModel
from obspy.imaging.cm import pqlx
from pycheron.util.format import parse_snclq
from pycheron.util.logger import Logger
from pycheron.util import *
from pycheron.metrics.psdMetric import psdMetric

from pycheron.db.sqllite_db import Database


def psdPlot(
    st,
    style="psd",
    f_name=None,
    showNoiseModel=True,
    showMaxMin=False,
    showMode=True,
    showMean=False,
    showMedian=False,
    showEnvelope=False,
    envelopeType="10_90",
    showSingle=False,
    singleType=None,
    pcolor=pqlx,
    ylo=-200,
    yhi=-50,
    timespan=0,
    evalresp=None,
    network=None,
    station=None,
    channel=None,
    location=None,
    session=None,
    database=None,
    logger=None,
):
    """
    Plot instrument corrected noise values of PSDs and PDFs

    :param st: (stream object or Database) - PSD or PDF as output from getNoise or getPDF or Database object
    :type st: obspy.core.stream.Stream or pycheron.db.sqllite_db.Database
    :param style: String to determine plotting style, options include psd or pdf
    :type style: str
    :param showNoiseModel: If True, display noise models on the plot
    :type showNoiseModel: bool
    :param showMaxMin: (boolean) - If True, display max and min on the plot
    :type showMaxMin: bool
    :param showMode: (boolean) - If True, display mode on the plot
    :type showMode: bool
    :param showMean: (boolean) - If True, display mean on the plot
    :type showMean: bool
    :param showMedian: (boolean) - If True, display median on the plot
    :type showMedian: bool
    :param showEnvelope: (boolean) - If true, percentile envelope will be plotted.
    :type showEnvelope: boo;
    :param envelopeType: (str) - envelope to be plotted (Default='10_90")
                "10_90" - 10th and 90th percentile envelope
                "05_95" - 5th annd 95th percentile envelope
    :type envelopeType: str
    :param showSingle: (bool) - if True, plot single percentile.
    :type showSingle: bool
    :param singleType: Single line percentile to be plotted. Options: "5" (5th percentile), "10" (10th percentile)
                       "90" ( 90th percentile), "95" (95th percentile).
    :type singleType: str
    :param pcolor: (colormap) - color map to utilize in plots
    :type pcolor: obspy.imaging.cm
    :param ylo: min dB to show on the y-axis scale
    :type ylo: int
    :param yhi: max dB to show on the y-axis scale
    :type yhi: int
    :param evalresp: evalresp directory or np.ndarray
    :type evalresp: str or numpy.ndarray
    :param timespan: (int) Time span for plots to be broken up into. 1 - 1 week, 2 - 2 weeks, 4 - 1 month, Any other int
                    will result in all data being processed.
    :type timespan: int
    :param network: If using database, network name
    :type network: str
    :param station: If using database, station name
    :type station: str
    :param channel: If using database, channel name
    :type channel: str
    :param location: If using database, location name
    :type location: str
    :param session: If using database, session name
    :type session: str
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger

    **Example**

    1. PSD Plot

    .. code-block:: python

        #PSD Plot

        import obspy
        from pycheron.psd.psdList import psdList
        from pycheron.psd.psdStatistics import psdStatistics
        from pycheron.plotting.psdPlot import psdPlot
        from pycheron.psd.getPDF import getPDF
        from pycheron.psd.noise.getNoise import getNoise
        from pycheron.psd.noise.noiseModel import noiseModel
        from obspy.imaging.cm import pqlx


        data = 'test/test_data/7a_cabn_bhe.884965.tar.mseed'

        #reading in stream - show PSD with Statistics and NHNM/NLNM
        st = obspy.read(data)


        psdPlot(st,style='psd', showNoiseModel=True, showMaxMin=True, showMode=True, showMean=True, showMedian=True, showEnvelope=True, envelopeType="10_90",pcolor=pqlx)

    .. image:: _static/psdPlot.png

    2. PDF Plot

    .. code-block:: python

        #PDF Plot

        import obspy
        from pycheron.psd.psdList import psdList
        from pycheron.psd.psdStatistics import psdStatistics
        from pycheron.plotting.psdPlot import psdPlot
        from pycheron.psd.getPDF import getPDF
        from pycheron.psd.noise.getNoise import getNoise
        from pycheron.psd.noise.noiseModel import noiseModel
        from obspy.imaging.cm import pqlx


        data = 'test/test_data/7a_cabn_bhe.884965.tar.mseed'

        #reading in stream - show PSD with Statistics and NHNM/NLNM
        st = obspy.read(data)


        psdPlot(st,style='pdf', showNoiseModel=True, showMaxMin=True, showMode=True, showMean=True, showMedian=True, showEnvelope=True, envelopeType="10_90",pcolor=pqlx)

    .. image:: _static/pdfPlot.png

    """

    if logger is None:
        logger = Logger(None)

    # Get basic sncl and time information from psd
    # If only one channel in file get the snclq information in the following manner
    # Register cmap
    plt.close()
    plt.style.use("default")
    plt.register_cmap(cmap=pqlx)
    if isinstance(st, Database):
        _psdPlot_from_database(
            st,
            timespan=timespan,
            style=style,
            f_name=f_name,
            showNoiseModel=showNoiseModel,
            showMaxMin=showMaxMin,
            showMode=showMode,
            showMean=showMean,
            showMedian=showMedian,
            showEnvelope=showEnvelope,
            envelopeType=envelopeType,
            showSingle=showSingle,
            singleType=singleType,
            pcolor=pcolor,
            network=network,
            station=station,
            channel=channel,
            location=location,
            session=session,
            logger=None,
        )
    else:
        try:
            network, station, channel, location, quality = parse_snclq(st[0].get_id())
        except:
            network, station, channel, location = parse_snclq(st[0].get_id())

        # If database object, check to see if we already have the psds calculated
        if database is not None:
            time = (st[0].stats.starttime.isoformat(), st[0].stats.endtime.isoformat())
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
            if tb.empty:
                PSD = psdList(st)
                psds_metrics_for_db = psdMetric(st)
                database.insert_metric(psds_metrics_for_db)

            else:
                tb_psds = tb.uncorrected_psds[tb.channel == channel]
                PSD = []
                for j in range(len(tb_psds)):
                    psd = database.extract_masks(tb_psds.iloc[j])
                    PSD.append(psd)
        else:
            PSD = psdList(st)

        if len(PSD) == 1:
            # for i in range(len(PSD[0])):
            #  #Count the number of PSD segments
            psdCount = len(PSD[0])

            # Get snclq information from psdList, then parse out sta, net, chan, loc, frequency.
            # Get starttime of first psd segment and endtime from last psd segment
            snclqSplit = (PSD[0][0][2][0]).split(".")
            network = snclqSplit[0]
            station = snclqSplit[1]
            channel = snclqSplit[3]
            starttime = PSD[0][0][3]
            endtime = PSD[0][-1][4]
            freq = np.asarray(PSD[0][0][0])

            # Generate basic stats as well as noiseMatrix
            stats = psdStatistics(
                PSD, evalresp=evalresp, logger=logger, database=database
            )
            noiseMatrix = stats[0]["noise_matrix_noise"]

            # Create period variable for determining plotting xlimits if not broadband B* channel name
            period = 1 / freq

            out = {
                "network": network,
                "station": station,
                "channel": channel,
                "location": location,
                "session": session,
                "metric_name": "psdPlot",
                "plot_style": style,
                "timespan": timespan,
            }

            # Default style is psd, but style may be changed to pdf (see below). Plot PSDS from noiseMatrix vs. period,
            # (x axis log), turn on grid and label axes. Have title show how many Corrected PSDs showing with snclq
            if style is "psd":
                for i in noiseMatrix:
                    plt.semilogx(period, i)
                plt.grid(True, which="both")
                plt.xlabel("Period (s)")
                plt.ylabel("Power (dB)")
                plt.title(
                    "{} Corrected, hourly PSDs from {}-{} for {}.{}.{}".format(
                        psdCount, starttime, endtime, network, station, channel
                    )
                )

            elif style is "pdf":
                # If style is pdf, flip PDFMatrix because using period, grab out PDF bins
                PDFMatrix = np.fliplr(stats[0]["pdf_matrix"])
                PDFBins = stats[0]["pdf_bins"]
                # Plot colormesh with X = rev(period), y=PDFBins, Z = PDFMatrix. Use pqlx color map and standard 0-30
                # scaling when showing McNamara PSD PDFs. Make x axis scaled logarithmically, label axes, turn on grid and
                # plot colorbar
                plt.pcolormesh(
                    period[::-1],
                    PDFBins,
                    PDFMatrix,
                    cmap=pcolor,
                    alpha=1.0,
                    vmin=0,
                    vmax=30,
                )
                plt.xscale("log")
                plt.ylabel("Amplitude [$m^2/s^4/Hz$] [dB]")
                plt.xlabel("Period (s)")
                plt.title(
                    "PDF plot of {} hourly PSDs from {}-{} for {}.{}.{}".format(
                        psdCount, starttime, endtime, network, station, channel
                    )
                )
                plt.grid(True, which="both")
                plt.colorbar()

            # Choose plotting limits based on channel names
            if channel.startswith("B"):
                plt.xlim(0.1, 100)
            else:
                if min(period) < 0.001:
                    plt.xlim(0.001, max(period))
                else:
                    plt.xlim(min(period), max(period))

            # Use dB limits in plot yaxis
            plt.ylim(ylo, yhi)

            # Add basic stat lines depending on above TRUE/FALSE call statements (e.g., showMode=TRUE)
            # Show NLNM/NHNM in grey if True; only show from 0.1 -100s
            if showNoiseModel:
                freq_sub = freq[freq <= 10]
                period_sub = 1 / freq[freq <= 10]
                nlnm, nhnm = noiseModel(freq_sub)
                plt.semilogx(period_sub, nlnm, linewidth=3, label="NLNM", color="grey")
                plt.semilogx(period_sub, nhnm, linewidth=3, label="NHNM", color="grey")
                plt.legend()
                update = {"type": "nlnm"}
                out.update(update)

            # Display max/Min if true as blue and red lines
            if showMaxMin:
                plt.semilogx(
                    period, stats[0]["max"], linewidth=3, label="Max", color="blue"
                )
                plt.semilogx(
                    period, stats[0]["min"], linewidth=3, label="Min", color="red"
                )
                plt.legend()
                update = {"type": "MinMax"}
                out.update(update)

            # Display mode if true in yellow
            if showMode:
                plt.semilogx(
                    period, stats[0]["mode"], linewidth=3, label="Mode", color="yellow"
                )
                plt.legend()
                update = {"type": "mode"}
                out.update(update)

            # Display mean if true in orange
            if showMean:
                plt.semilogx(
                    period, stats[0]["mean"], linewidth=3, label="Mean", color="orange"
                )
                plt.legend()
                update = {"type": "mean"}
                out.update(update)

            # Display median if true in green
            if showMedian:
                plt.semilogx(
                    period,
                    stats[0]["median"],
                    linewidth=3,
                    label="Median",
                    color="green",
                )
                plt.legend()
                update = {"type": "median"}
                out.update(update)

            # Display the percentile envelopes if true, 10 and 90 or 5 and 95
            if showEnvelope:
                if envelopeType == "10_90":
                    plt.semilogx(
                        period,
                        stats[0]["percent_10"],
                        linewidth=3,
                        label="10th Percentile",
                        color="black",
                    )
                    plt.semilogx(
                        period,
                        stats[0]["percent_90"],
                        linewidth=3,
                        label="90th Percentile",
                        color="black",
                    )
                    plt.legend()
                    update = {"type": "envelope_10_90"}
                    out.update(update)
                if envelopeType == "05_95":
                    plt.semilogx(
                        period,
                        stats[0]["percent_5"],
                        linewidth=3,
                        label="5th Percentile",
                        color="teal",
                    )
                    plt.semilogx(
                        period,
                        stats[0]["percent_95"],
                        linewidth=3,
                        label="95th Percentile",
                        color="teal",
                    )
                    plt.legend()
                    update = {"type": "envelope_05_95"}
                    out.update(update)

            # If true plot single percentile line
            if showSingle:
                plt.semilogx(
                    period,
                    stats[0]["percent_" + singleType],
                    linewidth=3,
                    label=singleType + "th Percentile",
                    color="c",
                )
                plt.legend()
                update = {"type": "percent_" + singleType}
                out.update(update)

            # Save figure into designated directory if f_name not none
            if f_name is None:
                f_name = "f_name"

            name = f_name + "_" + channel + ".png"
            plt.savefig(name)

            update = {"image_path": name}
            out.update(update)

            if isinstance(st, Database) or database is not None:
                database.insert_metric(out)

        # If more than one channel exists in the file, plotting and obtaining
        # snclq, stats, PSD information must account for all of the channels
        elif len(PSD) != 1:
            # # Create fig, ax handle for subplots defined by the number of PSDs, will allow user
            # # to see multiple channels on same plot
            fig, ax = plt.subplots(len(PSD), 1, figsize=(15, 15))

            # Cycle through the PSDs to get the appropriate snclq information from psdList,
            # then parse out sta,net,chan,loc,frequency. We expect these to be
            # #the same but in the case where a file might have multiple channel types we parse
            # them out individually in case the channels, locs, freq may be
            # different. Get starttime for first psd segment and endtime from last psd
            # segment for each different channel.

            for i in range(len(PSD)):
                try:
                    network, station, channel, location, quality = parse_snclq(
                        PSD[i][0][2][0]
                    )
                except Exception as e:
                    network, station, channel, location = parse_snclq(PSD[i][0][2][0])
                    print(
                        "Failed to parse snclq: {e} (PSD of size: {PSDsize})".format(
                            e=e, PSDsize=len(PSD)
                        )
                    )
                starttime = PSD[i][0][3]
                endtime = PSD[i][-1][4]
                freq = np.asarray(PSD[i][0][0])

                # Count the number of PSD segments

                psdCount = len(PSD[i])

                out = {
                    "network": network,
                    "station": station,
                    "channel": channel,
                    "location": location,
                    "session": session,
                    "metric_name": "psdPlot",
                    "plot_style": style,
                    "timespan": timespan,
                }

                # Generate basic stats as well as noiseMatrix.
                # Inside loop so that we can get the appropriate noise matrices for each channel
                stats = psdStatistics(
                    PSD, evalresp=evalresp, logger=logger, database=database
                )
                noiseMatrix = stats[i]["noise_matrix_noise"]

                # Create period variable for determining plotting xlimits if not broadband B*
                # channel name. On the off chance that the frequencies are different within the
                # file we will make a list of periods based on i.
                period = 1 / freq

                # Get noise and PDF based on psd
                f, n, PSD = getNoise(PSD, evalresp=None)

                # Default style is psd, but style may be changed to pdf (see below). Same as above: Plot PSDS from
                # noiseMatrix vs. period for each channel in PSD. (x axis log), turn on grid and label axes.
                # Have title show how many Corrected PSDs showing with snclq
                if style is "psd":
                    for j in noiseMatrix:
                        ax[i].semilogx(period, j)
                        ax[i].set_title(
                            "{} Corrected, hourly PSDs from {}-{} for {}.{}.{}".format(
                                psdCount, starttime, endtime, network, station, channel
                            )
                        )
                ax[i].grid(True, which="both")
                ax[i].set_xlabel("Period (s)")
                ax[i].set_ylabel("Power (dB)")

                if style is "pdf":
                    # If style is pdf, flip PDFMatrix because using period, grab out PDF bins
                    PDFMatrix = np.fliplr(stats[i]["pdf_matrix"])
                    PDFBins = stats[i]["pdf_bins"]
                    # Plot colormesh with X = rev(period), y=PDFBins, Z = PDFMatrix. Use pqlx color map and standard 0-30
                    # scaling when showing McNamara PSD PDFs. Make x axis scaled logarithmically, label axes, turn on grid and
                    # plot colorbar
                    im = ax[i].pcolormesh(
                        period[::-1],
                        PDFBins,
                        PDFMatrix,
                        cmap=pcolor,
                        alpha=1.0,
                        vmin=0,
                        vmax=30,
                    )
                    ax[i].set_xscale("log")
                    ax[i].set_ylabel("Amplitude [$m^2/s^4/Hz$] [dB]")
                    ax[i].set_xlabel("Period (s)")
                    ax[i].set_title(
                        "PDF plot of {} hourly PSDs from {}-{} for {}.{}.{}".format(
                            psdCount, starttime, endtime, network, station, channel
                        )
                    )
                    ax[i].grid(True, which="both")

                # Choose plotting limits based on channel names
                if channel.startswith("B"):
                    ax[i].set_xlim(0.1, 100)
                else:
                    if min(period) < 0.001:
                        ax[i].set_xlim(0.001, max(period))
                    else:
                        ax[i].set_xlim(min(period), max(period))

                # Use dB limits in plot yaxis
                ax[i].set_ylim(ylo, yhi)

                # Add basic stat lines depending on above TRUE/FALSE call statements (e.g., showMode=TRUE)
                #  Show NLNM/NHNM in grey if True; only show from 0.1 -100s
                if showNoiseModel:
                    freq_sub = freq[freq <= 10]
                    period_sub = 1 / freq[freq <= 10]
                    nlnm, nhnm = noiseModel(freq_sub)
                    ax[i].semilogx(
                        period_sub, nlnm, linewidth=3, label="NLNM", color="grey"
                    )
                    ax[i].semilogx(
                        period_sub, nhnm, linewidth=3, label="NHNM", color="grey"
                    )
                    ax[i].legend()
                    if i == 0:
                        update = {"type": "nlnm"}
                        out.update(update)

                # If True, display the min/max in blue, red
                if showMaxMin:
                    ax[i].semilogx(
                        period, stats[i]["max"], linewidth=3, label="Max", color="blue"
                    )
                    ax[i].semilogx(
                        period, stats[i]["min"], linewidth=3, label="Min", color="red"
                    )
                    ax[i].legend()
                    if i == 0:
                        update = {"type": "MinMax"}
                        out.update(update)

                # If True, display the mode in yellow
                if showMode:
                    ax[i].semilogx(
                        period,
                        stats[i]["mode"],
                        linewidth=3,
                        label="Mode",
                        color="yellow",
                    )
                    ax[i].legend()
                    if i == 0:
                        update = {"type": "mode"}
                        out.update(update)

                # If true, display the mean in orange
                if showMean:
                    ax[i].semilogx(
                        period,
                        stats[i]["mean"],
                        linewidth=3,
                        label="Mean",
                        color="orange",
                    )
                    ax[i].legend()
                    if i == 0:
                        update = {"type": "mean"}
                        out.update(update)

                # If true, display the median in green
                if showMedian:
                    ax[i].semilogx(
                        period,
                        stats[i]["median"],
                        linewidth=3,
                        label="Median",
                        color="green",
                    )
                    ax[i].legend()
                    if i == 0:
                        update = {"type": "median"}
                        out.update(update)

                # Display the percentile envelopes if true, 10 and 90 or 5 and 95
                if showEnvelope:
                    if envelopeType == "10_90":
                        ax[i].semilogx(
                            period,
                            stats[i]["percent_10"],
                            linewidth=3,
                            label="10th Percentile",
                            color="black",
                        )
                        ax[i].semilogx(
                            period,
                            stats[i]["percent_90"],
                            linewidth=3,
                            label="90th Percentile",
                            color="black",
                        )
                        ax[i].legend()
                        if i == 0:
                            update = {"type": "envelope_10_90"}
                            out.update(update)

                    if envelopeType == "05_95":
                        ax[i].semilogx(
                            period,
                            stats[i]["percent_5"],
                            linewidth=3,
                            label="5th Percentile",
                            color="teal",
                        )
                        ax[i].semilogx(
                            period,
                            stats[i]["percent_95"],
                            linewidth=3,
                            label="95th Percentile",
                            color="teal",
                        )
                        ax[i].legend()
                        if i == 0:
                            update = {"type": "envelope_05_95"}
                            out.update(update)

                # If true plot single percentile line
                if showSingle:
                    ax[i].semilogx(
                        period,
                        stats[i]["percent_" + singleType],
                        linewidth=3,
                        label=singleType + "th Percentile",
                        color="c",
                    )
                    ax[i].legend()
                    if i == 0:
                        update = {"type": "percent_" + singleType}
                        out.update(update)

            # Plot colorbar for pdfs
            if style is "pdf":
                cax, kw = mpl.colorbar.make_axes([axes for axes in ax.flat])
                plt.colorbar(im, cax=cax, **kw)
                cbar = plt.colorbar(im, cax=cax, **kw)
                cbar.set_label("Probability (%)")

            # Save figure into designated directory if f_name not none
            if f_name is None:
                f_name = "f_name"

            name = f_name + "_" + station + ".png"
            plt.savefig(name)

            update = {"image_path": name}
            out.update(update)

            if isinstance(st, Database) or database is not None:
                database.insert_metric(out)

            plt.close()


def get_psd_plot_data(
    database,
    style="psd",
    f_name=None,
    showNoiseModel=True,
    showMaxMin=False,
    showMode=False,
    showMean=False,
    evalresp=None,
    showMedian=True,
    showEnvelope=True,
    envelopeType="10_90",
    showSingle=False,
    singleType=None,
    ylo=-200,
    yhi=-50,
    pcolor=pqlx,
    timespan=1,
    network=None,
    station=None,
    channel=None,
    location=None,
    session=None,
    logger=None,
):
    """
    Internal function used by psdPlot if inout is Database
    """
    # Get basic sncl and time information from psd
    # If only one channel in file get the snclq information in the following manner
    psd_channels = _read_psds(
        database,
        timespan=timespan,
        network=network,
        station=station,
        channel=channel,
        location=location,
        session=session,
        logger=logger,
    )
    period = 0
    psdCount = 0
    starttime = 0
    endtime = 0
    freq = 0
    # If psd_channels is empty exit
    if not psd_channels:
        return

    PSD = []
    for ind in range(len(psd_channels)):
        PSD.append(psd_channels[ind][0][0][0])

    if len(PSD) == 1:
        PSD = psd_channels[0][0][0]

        psdCount = len(PSD[0])

        # Get snclq information from psdList, then parse out sta, net, chan, loc, frequency.
        # Get starttime of first psd segment and endtime from last psd segment
        snclqSplit = (PSD[0][0][2][0]).split(".")
        network = snclqSplit[0]
        station = snclqSplit[1]
        channel = snclqSplit[3]
        starttime = PSD[0][0][3]
        endtime = PSD[0][-1][4]
        freq = np.asarray(PSD[0][0][0])

        # Generate basic stats as well as noiseMatrix
        # stats = psdStatistics(PSD, evalresp=evalresp, logger=logger, database=database)
        # noiseMatrix = stats[0]["noise_matrix_noise"]

        # Create period variable for determining plotting xlimits if not broadband B* channel name
        stats = psdStatistics(PSD, evalresp=evalresp, logger=logger, database=database)
        noiseMatrix = stats[0]["noise_matrix_noise"]
        period = 1 / freq

    elif len(psd_channels) != 1:
        # # Create fig, ax handle for subplots defined by the number of PSDs, will allow user
        # # to see multiple channels on same plot
        fig, ax = plt.subplots(len(psd_channels), 1, figsize=(15, 15))

        # Cycle through the PSDs to get the appropriate snclq information from psdList,
        # then parse out sta,net,chan,loc,frequency. We expect these to be
        # #the same but in the case where a file might have multiple channel types we parse
        # them out individually in case the channels, locs, freq may be
        # different. Get starttime for first psd segment and endtime from last psd
        # segment for each different channel.

        for i in range(len(psd_channels)):
            try:
                network, station, channel, location, quality = parse_snclq(
                    PSD[i][0][2][0]
                )
            except:
                network, station, channel, location = parse_snclq(PSD[i][0][2][0])
            starttime = PSD[i][0][3]
            endtime = PSD[i][-1][4]
            freq = np.asarray(PSD[i][0][0])

            # Count the number of PSD segments

            psdCount = len(PSD[i])

            out = {
                "network": network,
                "station": station,
                "channel": channel,
                "location": location,
                "session": session,
                "metric_name": "psdPlot",
                "plot_style": style,
                "timespan": timespan,
            }

            # Generate basic stats as well as noiseMatrix.
            # Inside loop so that we can get the appropriate noise matrices for each channel
            stats = psdStatistics(
                PSD, evalresp=evalresp, logger=logger, database=database
            )
            noiseMatrix = stats[i]["noise_matrix_noise"]

            # Create period variable for determining plotting xlimits if not broadband B*
            # channel name. On the off chance that the frequencies are different within the
            # file we will make a list of periods based on i.

            period = 1 / freq

    return PSD, psdCount, period, style, stats, noiseMatrix, starttime, endtime, freq


def _psdPlot_from_database(
    database,
    style="psd",
    f_name=None,
    showNoiseModel=True,
    showMaxMin=False,
    showMode=False,
    showMean=False,
    evalresp=None,
    showMedian=True,
    showEnvelope=True,
    envelopeType="10_90",
    showSingle=False,
    singleType=None,
    ylo=-200,
    yhi=-50,
    pcolor=pqlx,
    timespan=1,
    network=None,
    station=None,
    channel=None,
    location=None,
    session=None,
    logger=None,
):
    """
    Internal function used by psdPlot if inout is Database
    """
    # Get basic sncl and time information from psd
    # If only one channel in file get the snclq information in the following manner

    psd_channels = _read_psds(
        database,
        timespan=timespan,
        network=network,
        station=station,
        channel=channel,
        location=location,
        session=session,
        logger=logger,
    )

    # If psd_channels is empty exit
    if not psd_channels:
        return

    PSD = []
    for ind in range(len(psd_channels)):
        PSD.append(psd_channels[ind][0][0][0])

    if len(PSD) == 1:
        PSD = psd_channels[0][0][0]

        psdCount = len(PSD[0])

        # Get snclq information from psdList, then parse out sta, net, chan, loc, frequency.
        # Get starttime of first psd segment and endtime from last psd segment
        snclqSplit = (PSD[0][0][2][0]).split(".")
        network = snclqSplit[0]
        station = snclqSplit[1]
        channel = snclqSplit[3]
        starttime = PSD[0][0][3]
        endtime = PSD[0][-1][4]
        freq = np.asarray(PSD[0][0][0])

        # Generate basic stats as well as noiseMatrix
        stats = psdStatistics(PSD, evalresp=evalresp, logger=logger, database=database)
        noiseMatrix = stats[0]["noise_matrix_noise"]

        # Create period variable for determining plotting xlimits if not broadband B* channel name
        period = 1 / freq

        out = {
            "network": network,
            "station": station,
            "channel": channel,
            "location": location,
            "session": session,
            "metric_name": "psdPlot",
            "plot_style": style,
            "timespan": timespan,
        }

        # Default style is psd, but style may be changed to pdf (see below). Plot PSDS from noiseMatrix vs. period,
        # (x axis log), turn on grid and label axes. Have title show how many Corrected PSDs showing with snclq
        if style is "psd":
            for i in noiseMatrix:
                plt.semilogx(period, i)
            plt.grid(True, which="both")
            plt.xlabel("Period (s)")
            plt.ylabel("Power (dB)")
            plt.title(
                "{} Corrected, hourly PSDs from {}-{} for {}.{}.{}".format(
                    psdCount, starttime, endtime, network, station, channel
                )
            )

        elif style is "pdf":
            # If style is pdf, flip PDFMatrix because using period, grab out PDF bins
            PDFMatrix = np.fliplr(stats[0]["pdf_matrix"])
            PDFBins = stats[0]["pdf_bins"]
            # Plot colormesh with X = rev(period), y=PDFBins, Z = PDFMatrix. Use pqlx color map and standard 0-30
            # scaling when showing McNamara PSD PDFs. Make x axis scaled logarithmically, label axes, turn on grid and
            # plot colorbar
            plt.pcolormesh(
                period[::-1],
                PDFBins,
                PDFMatrix,
                cmap=pcolor,
                alpha=1.0,
                vmin=0,
                vmax=30,
            )
            plt.xscale("log")
            plt.ylabel("Amplitude [$m^2/s^4/Hz$] [dB]")
            plt.xlabel("Period (s)")
            plt.title(
                "PDF plot of {} hourly PSDs from {}-{} for {}.{}.{}".format(
                    psdCount, starttime, endtime, network, station, channel
                )
            )
            plt.grid(True, which="both")
            plt.colorbar()

        # Choose plotting limits based on channel names
        if channel.startswith("B"):
            plt.xlim(0.1, 100)
        else:
            if min(period) < 0.001:
                plt.xlim(0.001, max(period))
            else:
                plt.xlim(min(period), max(period))

        # Use dB limits in plot yaxis
        plt.ylim(ylo, yhi)

        # Add basic stat lines depending on above TRUE/FALSE call statements (e.g., showMode=TRUE)
        # Show NLNM/NHNM in grey if True; only show from 0.1 -100s
        if showNoiseModel:
            freq_sub = freq[freq <= 10]
            period_sub = 1 / freq[freq <= 10]
            nlnm, nhnm = noiseModel(freq_sub)
            plt.semilogx(period_sub, nlnm, linewidth=3, label="NLNM", color="grey")
            plt.semilogx(period_sub, nhnm, linewidth=3, label="NHNM", color="grey")
            plt.legend()
            update = {"type": "nlnm"}
            out.update(update)

        # Display max/Min if true as blue and red lines
        if showMaxMin:
            plt.semilogx(
                period, stats[0]["max"], linewidth=3, label="Max", color="blue"
            )
            plt.semilogx(period, stats[0]["min"], linewidth=3, label="Min", color="red")
            plt.legend()
            update = {"type": "MinMax"}
            out.update(update)

        # Display mode if true in yellow
        if showMode:
            plt.semilogx(
                period, stats[0]["mode"], linewidth=3, label="Mode", color="yellow"
            )
            plt.legend()
            update = {"type": "mode"}
            out.update(update)

        # Display mean if true in orange
        if showMean:
            plt.semilogx(
                period, stats[0]["mean"], linewidth=3, label="Mean", color="orange"
            )
            plt.legend()
            update = {"type": "mean"}
            out.update(update)

        # Display median if true in green
        if showMedian:
            plt.semilogx(
                period, stats[0]["median"], linewidth=3, label="Median", color="green"
            )
            plt.legend()
            update = {"type": "median"}
            out.update(update)

        # Display the percentile envelopes if true, 10 and 90 or 5 and 95
        if showEnvelope:
            if envelopeType == "10_90":
                plt.semilogx(
                    period,
                    stats[0]["percent_10"],
                    linewidth=3,
                    label="10th Percentile",
                    color="black",
                )
                plt.semilogx(
                    period,
                    stats[0]["percent_90"],
                    linewidth=3,
                    label="90th Percentile",
                    color="black",
                )
                plt.legend()
                update = {"type": "envelope_10_90"}
                out.update(update)
            if envelopeType == "05_95":
                plt.semilogx(
                    period,
                    stats[0]["percent_5"],
                    linewidth=3,
                    label="5th Percentile",
                    color="teal",
                )
                plt.semilogx(
                    period,
                    stats[0]["percent_95"],
                    linewidth=3,
                    label="95th Percentile",
                    color="teal",
                )
                plt.legend()
                update = {"type": "envelope_05_95"}
                out.update(update)

        # If true plot single percentile line
        if showSingle:
            plt.semilogx(
                period,
                stats[0]["percent_" + singleType],
                linewidth=3,
                label=singleType + "th Percentile",
                color="c",
            )
            plt.legend()
            update = {"type": "percent_" + singleType}
            out.update(update)

        # Save figure into designated directory if f_name not none
        if f_name is None:
            f_name = "f_name"

        name = f_name + "_" + channel + ".png"
        plt.savefig(name)

        update = {"image_path": name}
        out.update(update)

        database.insert_metric(out)

    # If more than one channel exists in the file, plotting and obtaining
    # snclq, stats, PSD information must account for all of the channels
    elif len(psd_channels) != 1:
        # # Create fig, ax handle for subplots defined by the number of PSDs, will allow user
        # # to see multiple channels on same plot
        fig, ax = plt.subplots(len(psd_channels), 1, figsize=(15, 15))

        # Cycle through the PSDs to get the appropriate snclq information from psdList,
        # then parse out sta,net,chan,loc,frequency. We expect these to be
        # #the same but in the case where a file might have multiple channel types we parse
        # them out individually in case the channels, locs, freq may be
        # different. Get starttime for first psd segment and endtime from last psd
        # segment for each different channel.

        for i in range(len(psd_channels)):
            try:
                network, station, channel, location, quality = parse_snclq(
                    PSD[i][0][2][0]
                )
            except:
                network, station, channel, location = parse_snclq(PSD[i][0][2][0])
            starttime = PSD[i][0][3]
            endtime = PSD[i][-1][4]
            freq = np.asarray(PSD[i][0][0])

            # Count the number of PSD segments

            psdCount = len(PSD[i])

            out = {
                "network": network,
                "station": station,
                "channel": channel,
                "location": location,
                "session": session,
                "metric_name": "psdPlot",
                "plot_style": style,
                "timespan": timespan,
            }

            # Generate basic stats as well as noiseMatrix.
            # Inside loop so that we can get the appropriate noise matrices for each channel
            stats = psdStatistics(
                PSD, evalresp=evalresp, logger=logger, database=database
            )
            noiseMatrix = stats[i]["noise_matrix_noise"]

            # Create period variable for determining plotting xlimits if not broadband B*
            # channel name. On the off chance that the frequencies are different within the
            # file we will make a list of periods based on i.
            period = 1 / freq

            # Get noise and PDF based on psd
            f, n, PSD = getNoise(PSD, evalresp=None)

            # Default style is psd, but style may be changed to pdf (see below). Same as above: Plot PSDS from
            # noiseMatrix vs. period for each channel in PSD. (x axis log), turn on grid and label axes.
            # Have title show how many Corrected PSDs showing with snclq
            if style is "psd":
                for j in noiseMatrix:
                    ax[i].semilogx(period, j)
                    ax[i].set_title(
                        "{} Corrected, hourly PSDs from {}-{} for {}.{}.{}".format(
                            psdCount, starttime, endtime, network, station, channel
                        )
                    )
            ax[i].grid(True, which="both")
            ax[i].set_xlabel("Period (s)")
            ax[i].set_ylabel("Power (dB)")

            if style is "pdf":
                # If style is pdf, flip PDFMatrix because using period, grab out PDF bins
                PDFMatrix = np.fliplr(stats[i]["pdf_matrix"])
                PDFBins = stats[i]["pdf_bins"]
                # Plot colormesh with X = rev(period), y=PDFBins, Z = PDFMatrix. Use pqlx color map and standard 0-30
                # scaling when showing McNamara PSD PDFs. Make x axis scaled logarithmically, label axes, turn on grid and
                # plot colorbar
                im = ax[i].pcolormesh(
                    period[::-1],
                    PDFBins,
                    PDFMatrix,
                    cmap=pcolor,
                    alpha=1.0,
                    vmin=0,
                    vmax=30,
                )
                ax[i].set_xscale("log")
                ax[i].set_ylabel("Amplitude [$m^2/s^4/Hz$] [dB]")
                ax[i].set_xlabel("Period (s)")
                ax[i].set_title(
                    "PDF plot of {} hourly PSDs from {}-{} for {}.{}.{}".format(
                        psdCount, starttime, endtime, network, station, channel
                    )
                )
                ax[i].grid(True, which="both")

            # Choose plotting limits based on channel names
            if channel.startswith("B"):
                ax[i].set_xlim(0.1, 100)
            else:
                if min(period) < 0.001:
                    ax[i].set_xlim(0.001, max(period))
                else:
                    ax[i].set_xlim(min(period), max(period))

            # Use dB limits in plot yaxis
            ax[i].set_ylim(ylo, yhi)

            # Add basic stat lines depending on above TRUE/FALSE call statements (e.g., showMode=TRUE)
            #  Show NLNM/NHNM in grey if True; only show from 0.1 -100s
            if showNoiseModel:
                freq_sub = freq[freq <= 10]
                period_sub = 1 / freq[freq <= 10]
                nlnm, nhnm = noiseModel(freq_sub)
                ax[i].semilogx(
                    period_sub, nlnm, linewidth=3, label="NLNM", color="grey"
                )
                ax[i].semilogx(
                    period_sub, nhnm, linewidth=3, label="NHNM", color="grey"
                )
                ax[i].legend()
                if i == 0:
                    update = {"type": "nlnm"}
                    out.update(update)

            # If True, display the min/max in blue, red
            if showMaxMin:
                ax[i].semilogx(
                    period, stats[i]["max"], linewidth=3, label="Max", color="blue"
                )
                ax[i].semilogx(
                    period, stats[i]["min"], linewidth=3, label="Min", color="red"
                )
                ax[i].legend()
                if i == 0:
                    update = {"type": "MinMax"}
                    out.update(update)

            # If True, display the mode in yellow
            if showMode:
                ax[i].semilogx(
                    period, stats[i]["mode"], linewidth=3, label="Mode", color="yellow"
                )
                ax[i].legend()
                if i == 0:
                    update = {"type": "mode"}
                    out.update(update)

            # If true, display the mean in orange
            if showMean:
                ax[i].semilogx(
                    period, stats[i]["mean"], linewidth=3, label="Mean", color="orange"
                )
                ax[i].legend()
                if i == 0:
                    update = {"type": "mean"}
                    out.update(update)

            # If true, display the median in green
            if showMedian:
                ax[i].semilogx(
                    period,
                    stats[i]["median"],
                    linewidth=3,
                    label="Median",
                    color="green",
                )
                ax[i].legend()
                if i == 0:
                    update = {"type": "median"}
                    out.update(update)

            # Display the percentile envelopes if true, 10 and 90 or 5 and 95
            if showEnvelope:
                if envelopeType == "10_90":
                    ax[i].semilogx(
                        period,
                        stats[i]["percent_10"],
                        linewidth=3,
                        label="10th Percentile",
                        color="black",
                    )
                    ax[i].semilogx(
                        period,
                        stats[i]["percent_90"],
                        linewidth=3,
                        label="90th Percentile",
                        color="black",
                    )
                    ax[i].legend()
                    if i == 0:
                        update = {"type": "envelope_10_90"}
                        out.update(update)

                if envelopeType == "05_95":
                    ax[i].semilogx(
                        period,
                        stats[i]["percent_5"],
                        linewidth=3,
                        label="5th Percentile",
                        color="teal",
                    )
                    ax[i].semilogx(
                        period,
                        stats[i]["percent_95"],
                        linewidth=3,
                        label="95th Percentile",
                        color="teal",
                    )
                    ax[i].legend()
                    if i == 0:
                        update = {"type": "envelope_05_95"}
                        out.update(update)

            # If true plot single percentile line
            if showSingle:
                ax[i].semilogx(
                    period,
                    stats[i]["percent_" + singleType],
                    linewidth=3,
                    label=singleType + "th Percentile",
                    color="c",
                )
                ax[i].legend()
                if i == 0:
                    update = {"type": "percent_" + singleType}
                    out.update(update)

        # Plot colorbar for pdfs
        if style is "pdf":
            cax, kw = mpl.colorbar.make_axes([axes for axes in ax.flat])
            plt.colorbar(im, cax=cax, **kw)
            cbar = plt.colorbar(im, cax=cax, **kw)
            cbar.set_label("Probability (%)")

        # Save figure into designated directory if f_name not none
        if f_name is None:
            f_name = "f_name"

        name = f_name + "_" + station + ".png"
        plt.savefig(name)

        update = {"image_path": name}
        out.update(update)

        database.insert_metric(out)

        plt.close()


def _read_psds(
    database,
    timespan=0,
    network=None,
    station=None,
    channel=None,
    location=None,
    session=None,
    logger=None,
):
    """
    Internal function to read psd from db
    """
    if logger is None:
        logger = Logger(None)

    psds = []
    chans = []
    d = []

    tb = database.get_metric(
        "psdMetric", network, station, channel, location, session=session
    )

    if tb.empty:
        logger.error(
            "psdPlot Error: No psdMetric found for specific network in Database, please run psdMetric first"
        )
        return

    else:
        tb_psds = tb.uncorrected_psds[tb.station == station]
        for j in range(len(tb_psds)):
            psd = database.extract_masks(tb_psds.iloc[j])
            try:
                snclq = psd[0][2][0]
            except IndexError:
                continue
            channel = snclq.split(".")[3]
            chans.append(channel)
            psds.append(psd)
            d.append(psd[0][3])
    chans = np.unique(chans)
    d = sorted(np.unique(d))
    # Setting date ranges
    all_dates = []
    dates = []
    # one week of data
    if timespan == 1:
        i = 0
        while len(dates) <= 6 and i < len(d):
            dates.append(d[i])
            i += 1
            if len(dates) == 7 or i == len(d):
                all_dates.append(dates)
                dates = []

    # 2 weeks of data
    elif timespan == 2:
        i = 0
        while len(dates) <= 13 and i < len(d):
            dates.append(d[i])
            i += 1
            if len(dates) == 14 or i == len(d):
                all_dates.append(dates)
                dates = []

    # one month of data
    elif timespan == 4:
        i = 0
        while len(dates) <= 27 and i < len(d):
            dates.append(d[i])
            i += 1
            if len(dates) == 28 or i == len(d):
                all_dates.append(dates)
                dates = []
    # plot all
    else:
        all_dates.append(d)

    all = []
    # loop through channels
    for chan in chans:
        # List that holds a list of date ranges, containing a list of channels for each date range. ex.
        # [7-day date range1:[ chan1: [day 1: [psds], day 2: [psds]], chan2: [etc.]],
        # 7-day date range1:[ chan1: [day 8: [psds], day 9: [psds], chan2: [etc.]]]

        all_subchan = []
        # loop thru the date ranges
        for date_range in all_dates:
            date_chunk = []
            chan_e = []
            chan_n = []
            chan_z = []
            for day in date_range:
                for psd in psds:
                    try:
                        snclq = psd[0][2][0]
                    except IndexError:
                        print(
                            "dailyPDFPlot Warning: Error in PSD file. Skipping PSD..."
                        )
                        continue
                    channel = snclq.split(".")[3]
                    start = psd[0][3]
                    # if first letter matches
                    if channel.startswith(chan) and start == day:
                        if channel.endswith("E") or channel.endswith("1"):
                            chan_e.append(psd)
                        elif channel.endswith("N") or channel.endswith("2"):
                            chan_n.append(psd)
                        else:
                            chan_z.append(psd)
                # making sure channel is not empty and not missing sub channel (ei E, N, Z)
            if chan_e:
                date_chunk.append(chan_e)
            if chan_n:
                date_chunk.append(chan_n)
            if chan_z:
                date_chunk.append(chan_z)

            all_subchan.append(date_chunk)
        all.append(all_subchan)
    return all
