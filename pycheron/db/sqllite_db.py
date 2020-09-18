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

import io
import codecs
import csv
import json
import os
import os.path
import sqlite3
import time
from sqlite3 import OperationalError, DatabaseError

import numpy as np
import pandas as pd
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client

from pycheron.util.format import *
from pycheron.util.getLatLon import get_latlon

__all__ = ["Database", "UnicodeWriter"]


class Database:
    """
    This is a custom database class which interfaces with callPycheronMetric as well as the other metric functions in
    the library. The database creates (if not already created) an SQLlite3 database and connect to it.

    :param db_name: Name of database
    :type db_name: str
    :param session_name: Session name. Allows for multiple sessions for running the same data
    :type session_name: str
    :param overwrite: If `True` any data that is the same will be overwritten
    :type overwrite: bool

    To create a database:

    .. code-block:: python

        >>> db = Database("test.db", session_name="session1", overwrite=True)


    """

    def __init__(self, db_name, session_name=None, overwrite=True):
        self._name = db_name
        self._session_name = session_name
        self._overwrite = overwrite
        try:
            self._conn = sqlite3.connect(db_name)

            db = self._conn.cursor()

            # make main table
            db.execute(
                """create table pycheron (
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               start_time varchar,
                               end_time varchar,
                               metric VARCHAR,
                               latitude varchar,
                               longitude varchar
                               )"""
            )

            # make basicStatsMetric
            db.execute(
                """create table basicStatsMetric (
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               metric_name varchar,
                               snclq varchar,
                               number_of_samples varchar,
                               start_time varchar,
                               end_time varchar,
                               minimum varchar,
                               maximum varchar,
                               median varchar,
                               mean varchar,
                               variance varchar,
                               standard_deviation varchar,
                               rms varchar,
                               min_mask varchar,
                               max_mask varchar,
                               median_mask varchar,
                               mean_mask varchar,
                               variance_mask varchar,
                               std_mask varchar,
                               rms_mask varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))"""
            )

            # make correlationMetric table
            db.execute(
                """create table correlationMetric (
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               correlation_coefficient varchar,
                               p_value varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))"""
            )

            # make crossCorrMetric table
            db.execute(
                """create table crossCorrMetric (
                                created varchar primary key,
                                network varchar,
                                station varchar,
                                channel varchar,
                                location varchar,
                                session varchar,
                                snclq varchar,
                                start_time varchar,
                                end_time varchar,
                                metric_name varchar,
                                peak_correlation varchar,
                                peak_lag varchar,
                                FOREIGN KEY (metric_name) REFERENCES pycheron(metric))
            """
            )

            # make DCOffSetMetricTimes table
            db.execute(
                """create table DCOffsetTimesMetric (
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               dc_offset_times varchar,
                               masks varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))"""
            )

            # make dedChanADFMetric table
            db.execute(
                """create table deadChannel(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               dead_channel_times varchar,
                               masks varchar,
                               is_dead_channel varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))"""
            )

            # make gapMetricSummary table
            db.execute(
                """create table gapMetricStation(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               session varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               station_completeness varchar,
                               total_gaps varchar,
                               total_overlaps varchar,
                               gap_masks varchar,
                               overlap_masks varchar,
                               combined_masks varchar,
                               maximum_overlap varchar,
                               maximum_gap varchar,
                               channel_percent_available varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station))"""
            )

            # make gapMetricChannel table
            db.execute(
                """create table gapMetricChannel(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               metric_name varchar,
                               duration varchar,
                               gap_overlap_start_time varchar,
                               gap_overlap_end_time varchar,
                               type varchar,
                               FOREIGN KEY (network) REFERENCES gapMetricStation(network),
                               FOREIGN KEY (station) REFERENCES gapMetricStation(station),
                               FOREIGN KEY (channel) REFERENCES gapMetricStation(channel))"""
            )

            # make psdMetric table
            db.execute(
                """create table psdMetric(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               noise2_mask varchar,
                               dead_channel_exponent varchar,
                               percent_above_nhnm varchar,
                               dead_channel_linear varchar,
                               uncorrected_psds varchar,
                               percent_below_nlnm varchar,
                               noise1_mask varchar,
                               bad_resp_mask varchar,
                               pdfs varchar,
                               hi_amp_mask varchar,
                               dc_mask varchar,
                               low_amp_mask varchar,
                               corrected_psds varchar,
                               dead_channel_gsn varchar,
                               dead_channel_exponent_hourly varchar,
                               dead_channel_linear_hourly varchar,
                               dead_channel_gsn_hourly varchar,
                               dead_chan_exp_hourly_masks varchar,
                               dead_chan_lin_hourly_masks varchar,
                               dead_chan_gsn_hourly_masks varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))"""
            )

            # make repeatedAmplitudeMetric table
            db.execute(
                """create table repeatedAmplitudeMetric(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               count varchar,
                               mask varchar,
                               rep_amp varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))"""
            )

            # make snrMetric table
            db.execute(
                """create table snrMetric(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               snr varchar,
                               algorithm varchar,
                               masks varchar,
                               windowSecs varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))"""
            )

            # make sohMetricActivityFlags table
            db.execute(
                """ create table sohMetricActivityFlags (
                           created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               session varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               calibration_signal varchar,
                               event_begin varchar,
                               event_end varchar,
                               event_in_progress varchar,
                               negative_leap varchar,
                               positive_leap varchar,
                               time_correction_applied varchar,                    
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station))

               """
            )

            # make sohMetricIOClockFlags table
            db.execute(
                """ create table sohMetricIOClockFlags (
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               session varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               clock_locked varchar,
                               end_time_series varchar,
                               long_record_read varchar,
                               short_record_read varchar,
                               start_time_series varchar,
                               station_volume varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station))

               """
            )
            # make sohMetricDataQualityFlags table
            db.execute(
                """create table sohMetricDataQualityFlags (
                           created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               session varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               amplifier_saturation varchar,
                               digital_filter_charging varchar,
                               digitizer_clipping varchar,
                               glitches varchar,
                               missing_padded_data varchar,
                               spikes varchar,
                               suspect_time_tag varchar,
                               telemetry_sync_error varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station))

               """
            )

            # make spikesMetric table
            db.execute(
                """create table spikesMetric(
               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               non_adjacent_spikes varchar,
                               mask varchar,
                               spike_times varchar, 
                               total_spike_count varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))

               """
            )

            # make staltaMetric table
            db.execute(
                """create table staltaMetric(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               event_time varchar,
                               max_stalta varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))

               """
            )

            # make transferFucntionMetric table
            db.execute(
                """create table transferFunctionMetric(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               gain_ratio varchar,
                               phase_diff varchar,
                               ms_coherence varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))

               """
            )

            # make dailyPDFPlot table
            db.execute(
                """create table dailyPdfPlot(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               channel varchar,
                               location varchar,
                               session varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               metric_name varchar,
                               noise_masks varchar,
                               microseism_masks varchar,
                               banded_masks varchar,
                               pdf_grid_image_path varchar,
                               pdf_line_image_path varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))

               """
            )

            # make networkNoiseModel table
            db.execute(
                """create table networkNoiseModel(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               session varchar,
                               metric_name varchar,
                               enModel varchar,
                               zModel varchar,
                               stations_en varchar,
                               stations_z varchar,
                               psdStats varchar,
                               network_noise_image_path varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network))                 
            """
            )

            # make stationNoiseModel table
            db.execute(
                """create table stationNoiseModel(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               session varchar,
                               metric_name varchar,
                               type varchar,
                               station_noise_image_path varchar,
                               e_station_noiseModel varchar,
                               n_station_noiseModel varchar,
                               z_station_noiseModel varchar,
                               e_station_noiseModel_low varchar,
                               n_station_noiseModel_low varchar,
                               z_station_noiseModel_low varchar, 
                               e_station_noiseModel_high varchar,
                               n_station_noiseModel_high varchar,
                               z_station_noiseModel_high varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station))
            """
            )

            # make psdPlot table
            db.execute(
                """create table psdPlot(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               channel varchar,
                               session varchar,
                               metric_name varchar,
                               plot_style varchar,
                               type varchar,
                               timespan varchar,
                               image_path varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))"""
            )

            # make calibrationMetric table
            db.execute(
                """create table calibrationMetric(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               channel varchar,
                               session varchar,
                               metric_name varchar,
                               start_time varchar,
                               end_time varchar,
                               snclq varchar,
                               detections varchar,
                               num_cals_detected varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))"""
            )

            # make stationRankingPlot table
            db.execute(
                """create table stationRankingPlot(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               channel varchar,
                               session varchar,
                               metric_name varchar,
                               rank_by varchar,
                               image_path varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel)
            )"""
            )

            # make psdStatistics table
            db.execute(
                """create table psdStatistics(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               channel varchar,
                               session varchar,
                               metric_name varchar,
                               noise_matrix_frequency varchar,
                               noise_matrix_noise varchar,
                               pdf_matrix varchar,
                               pdf_bins varchar,
                               max varchar,
                               min varchar,
                               mean varchar,
                               median varchar,
                               mode varchar,
                               nlnm varchar,
                               nhnm varchar,
                               percent_above_nhnm varchar,
                               percent_below_nlnm varchar,
                               percent_95 varchar,
                               percent_90 varchar,
                               percent_10 varchar,
                               percent_5 varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar, 
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel)
            )"""
            )

        except (OperationalError, ValueError):
            # if the database already exists, connect
            print(
                "{name} already exists. Connecting to {name}...".format(name=self._name)
            )
            self._conn = self.connect(db_name)

    def set_overwrite(self, overwrite):
        """Sets the overwrite parameters for database

        :param overwrite: If True, entries with the same session name, SNCLQ, and time will be over written.
        :type overwrite: bool

        """
        self._overwrite = overwrite

    def set_session_name(self, session_name):
        """Sets the session name. Purpose of this is to separate different runs with the same data

        :param session_name: session name
        :type session_name: str

        """
        self._session_name = session_name

    def get_connection(self):
        """Gets the database connection.

        :return: database connection

        """
        return self._conn

    def connect(self, db_name):
        """Connects to a specific database.

        :param db_name: database name
        :type db_name: str

        :return: database connection

        """
        self._name = db_name
        conn = sqlite3.connect(db_name)
        return conn

    def unique_with_nones(self, nparr):
        """Grabs unique elements from a numpy array, including None types
        :param nparr: numpy array of values
        :type nparr: numpy.ndarray

        :return: np.ndarray
        """
        if None in nparr:
            no_nones = nparr[nparr != np.array(None)]
            unique_no_nones = np.unique(no_nones)
            unique_final = np.append(unique_no_nones, None)
        else:
            unique_final = np.unique(nparr)
        return unique_final

    def networks(self):
        """Lists all of the unique networks in the database.

        :return: network names
        :rtype: list

        """
        db = self._conn.cursor()
        db.execute("""select ({cn}) from {tn}""".format(tn="pycheron", cn="network"))
        return list(self.unique_with_nones(np.array(self._extract(db.fetchall()))))

    def stations(self):
        """Lists all of the unique stations in the database.

        :return: station names
        :rtype: list

        """
        db = self._conn.cursor()
        db.execute("""select ({cn}) from {tn}""".format(tn="pycheron", cn="station"))
        return list(self.unique_with_nones(np.array(self._extract(db.fetchall()))))

    def channels(self):
        """Lists all of the unique channels in the database.

        :return: channel names
        :rtype: list

        """
        db = self._conn.cursor()
        db.execute("""select ({cn}) from {tn}""".format(tn="pycheron", cn="channel"))
        return list(self.unique_with_nones(np.array(self._extract(db.fetchall()))))

    def metrics(self):
        """Lists all of the unique metrics in the database.

        :return: metric names
        :rtype: list

        """
        db = self._conn.cursor()
        db.execute("""select ({cn}) from {tn}""".format(tn="pycheron", cn="metric"))
        return list(self.unique_with_nones(np.array(self._extract(db.fetchall()))))

    def session(self):
        """Lists all of the sessions in the database

        :return: session names
        :rtype: list

        """
        db = self._conn.cursor()
        db.execute("""select ({cn}) from {tn}""".format(tn="pycheron", cn="session"))
        return list(self.unique_with_nones(np.array(self._extract(db.fetchall()))))

    def view(
        self,
        metric_name=None,
        network=None,
        station=None,
        channel=None,
        location=None,
        time=None,
        session=None,
    ):
        """View in the main directory table (pycheron) all of the databases contents. Sort by combinations of metric name, network, station, channel, locations, time, and session name. If no options are given, it will return everything in the pycheron table.

        :param metric_name: metric name
        :type metric_name: str
        :param network: network name
        :type network: str
        :param station: station name
        :type station: str
        :param channel: channel name
        :type channel: str
        :param location: location
        :type location: str
        :param time: Tuple of starttime and endtime strings. ex. `time = ("2011-01-01T00:00:00","2011-01-01T23:55:13.190000")`
        :type time: tuple
        :param session: session name
        :type session: str

        :return: Pandas DataFrame
        :rtype: `pandas.DataFrame`
        """
        # if no location, change to --, otherwise it will be unsearchable/filterable
        if location == None or location == "":
            location = "--"

        # dictionary of parameters
        d = {
            "metric": metric_name,
            "network": network,
            "station": station,
            "location": location,
            "channel": channel,
            "session": session,
            "time": time,
        }
        # gets everything from pycheron table
        try:
            tb = pd.read_sql_query(
                """select {cn} from {tn}""".format(tn="pycheron", cn="*"), self._conn
            )
        # if gets error, try again
        except DatabaseError:
            tb = pd.read_sql_query(
                """select {cn} from {tn}""".format(tn="pycheron", cn="*"), self._conn
            )

        # loop through parameters dictionary
        for key, value in d.items():
            # if there is no value, skip. If the key is time, skip
            if value is not None and key != "time":
                # grabs the location in the database table where it matches the input param
                tb = tb.loc[tb[key] == value]
            # if the key is time, and not none
            elif value is not None and key == "time":
                # split the tuple into start and end time
                st = value[0]
                et = value[1]

                # edge case for crossCorrMetric, still todo
                if metric_name == "crossCorrMetric":
                    pass
                    # tb = tb.loc[tb["start_time"] >= str(UTCDateTime(st))]  # todo
                    # tb = tb.loc[tb["end_time"] <= str(UTCDateTime(et))]  # todo
                else:
                    tb = tb.loc[tb["start_time"] >= st]
                    tb = tb.loc[tb["end_time"] <= et]
        # if no results are returned.
        if tb.empty:
            print("No results found")
            return tb
        else:
            return tb

    def get_metric(
        self,
        metric_name,
        network=None,
        station=None,
        channel=None,
        location=None,
        time=None,
        session=None,
    ):
        """Gets data for specific metric. Searchable by network, station, channel, location, time and session name

        :param metric_name: metric name
        :type metric_name: str
        :param network: network name
        :type network: str
        :param station: station name
        :type station: str
        :param channel: channel name
        :type channel: str
        :param location: location
        :type location: str
        :param time: Tuple of starttime and endtime. ex. `time = ("2011-01-01T00:00:00","2011-01-01T23:55:13.190000")`
        :type time: tuple
        :param session: session name
        :type session: str

        :return: Pandas DataFrame
        :rtype: `pandas.DataFrame`

        """
        # connects to database
        db = self._conn.cursor()

        # if there is no location, set to --, otherwise it is unsearchable. networkNoiseModel does not have a location
        # so skip
        if (location == None or location == "") and metric_name != "networkNoiseModel":
            location = "--"

        # param dictionary
        d = {
            "network": network,
            "station": station,
            "location": location,
            "channel": channel,
            "time": time,
            "session": session,
        }
        # gets everything from specific metric table
        tb = pd.read_sql_query(
            """select {cn} from {tn}""".format(tn=metric_name, cn="*"), self._conn
        )

        # loops thru input params
        for key, value in d.items():
            # if value is not set, or key is time, skip
            if value is not None and key != "time":
                # matching input param to database values
                tb = tb.loc[tb[key] == value]
            # if key is time, split into start and end time, then filter.
            elif value is not None and key == "time":
                st = value[0]
                et = value[1]
                tb = tb.loc[tb["start_time"] >= st]
                tb = tb.loc[tb["end_time"] <= et]

        # if no results are returned
        if tb.empty:
            print("No results found")
            return tb
        else:
            return tb

    def query(self, sql):
        """For advanced users, a direct SQL command or statement can be used.

        :param sql: SQL formated statement
        :type sql: str

        :return: Pandas DataFrame
        :rtype: `pandas.DataFrame`
        """
        tb = pd.read_sql_query(sql, self._conn)
        return tb

    def insert_metric(self, metric):
        """Inserts metric into database

        :param metric: metric results
        :type metric: data

        """
        # get IRIS client
        client = Client("IRIS")
        # get connection
        conn = self._conn
        # connect to database
        db = conn.cursor()
        # get session name
        session = self._session_name

        # try/except sequence that deals with metrics whose names are inside lists or lists of lists
        try:
            metric_name = metric["metric_name"]
        except TypeError:
            try:
                metric_name = metric[0]["metric_name"]
            except TypeError:
                metric_name = metric[0][0]["metric_name"]

        # deals with metrics that have two SNCLQs
        if (
            metric_name == "correlationMetric"
            or metric_name == "crossCorrMetric"
            or metric_name == "transferFunctionMetric"
        ):
            # getting snclqs for the data
            snclq_combined = metric["snclq"]
            # parses snclq which are formatted: SNCLQ_1:SNCLQ_2
            snclq1, snclq2 = parse_snclq(snclq_combined, comb=True)

            # parses SNCLQs into individual parts.
            network, station, channel, location, quality = parse_snclq(snclq1)
            network2, station2, channel2, location2, quality2 = parse_snclq(snclq2)

            # if location is not, change to -- otherwise it is unsortable
            if location == "" and location2 == "":
                location = "--"
                location2 = "--"

            # Create the timestamp for the metric. This will be used as the primary key.
            timestamp = UTCDateTime(time.time()).isoformat()

            # getting the starttime and endtimes for the metric
            if metric_name == "crossCorrMetric":
                st = metric["snclq1_starttime"] + ":" + metric["snclq2_starttime"]
                et = metric["snclq1_endtime"] + ":" + metric["snclq2_endtime"]

                lat, lon = get_latlon(
                    client,
                    network,
                    station,
                    metric["snclq1_starttime"],
                    metric["snclq1_endtime"],
                )

                lat2, lon2 = get_latlon(
                    client,
                    network2,
                    station2,
                    metric["snclq2_starttime"],
                    metric["snclq2_endtime"],
                )

                lat = "{lat}:{lat2}".format(lat=lat, lat2=lat2)
                lon = "{lon}:{lon2}".format(lon=lon, lon2=lon2)
            else:
                st = metric["start_time"]
                et = metric["end_time"]

                lat, lon = get_latlon(
                    client,
                    network,
                    station,
                    st,
                    et,
                )

            # formatting the params to be inserted.
            net = "{n}:{n2}".format(n=network, n2=network2)
            sta = "{s}:{s2}".format(s=station, s2=station2)
            chan = "{c}:{c2}".format(c=channel, c2=channel2)
            loc = "{l}:{l2}".format(l=location, l2=location2)

            # checks current db to see if metric is already in db
            current_db = self.view(
                metric_name=metric_name,
                network=net,
                station=sta,
                channel=chan,
                location=loc,
                session=session,
                time=(st, et),
            )

            # if it is not currently in db
            if current_db.empty:
                names = """INSERT INTO pycheron( network, station, channel, location, session, metric, created, start_time, end_time, latitude, longitude )
                 VALUES ('{n}','{s}', '{c}', '{l}', '{sn}', '{metric}', '{created}', '{st}', '{et}', '{lat}', '{lon}')""".format(
                    n=net,
                    s=sta,
                    c=chan,
                    l=loc,
                    sn=session,
                    metric=metric_name,
                    created=timestamp,
                    st=st,
                    et=et,
                    lat=lat,
                    lon=lon,
                )
                db.execute(names)
                self._insert_metric_wrapper(conn, metric, timestamp)

            # if it is in the database, check the overwrite param, if true, overwrite.
            elif self._overwrite:
                new_vals = [
                    timestamp,
                    net,
                    sta,
                    chan,
                    loc,
                    session,
                    st,
                    et,
                    metric_name,
                    lat,
                    lon,
                ]
                idx = list(current_db.created.to_dict().keys())[0]
                old_created = current_db.created.to_dict()[idx]
                cn = list(current_db.keys())

                for i in range(len(new_vals)):
                    names = """UPDATE pycheron SET {cn} = '{new}' where created = '{c}'""".format(
                        cn=cn[i], new=new_vals[i], c=old_created
                    )
                    db.execute(names)

                self._insert_metric_wrapper(
                    conn, metric, timestamp, overwrite=old_created
                )

            else:
                print("Entry already exists, turn on overwrite or change session name")

        # deals with gapMetric
        elif "gapMetric" in metric_name:
            for i in range(len(metric[0])):
                metricStation = metric[0][i]
                timestamp = UTCDateTime(time.time()).isoformat()
                network = metricStation["network"]
                station = metricStation["station"]
                location = metricStation["location"]
                st = metricStation["start_time"].isoformat()
                et = metricStation["end_time"].isoformat()
                metric_name = metricStation["metric_name"]
                lat, lon = get_latlon(client, network, station, st, et)
                channel = None

                if location is None or location == "":
                    location = "--"

                current_db = self.view(
                    metric_name,
                    network,
                    station,
                    channel,
                    location,
                    session=session,
                    time=(st, et),
                )

                # if not currently in db
                if current_db.empty:
                    names = """INSERT INTO pycheron( network, station, channel, location, session, metric, created,
                     start_time, end_time, latitude, longitude ) VALUES ('{n}','{s}', '{c}', '{l}', '{sn}', '{metric}', '{created}', '{st}', '{et}', '{lat}', '{lon}')""".format(
                        n=network,
                        s=station,
                        c=channel,
                        metric=metric_name,
                        created=timestamp,
                        l=location,
                        sn=self._session_name,
                        st=st,
                        et=et,
                        lat=lat,
                        lon=lon,
                    )
                    db.execute(names)

                    self._insert_metric_wrapper(conn, metricStation, timestamp)

                elif self._overwrite:
                    new_vals = [
                        timestamp,
                        network,
                        station,
                        channel,
                        location,
                        session,
                        st,
                        et,
                        metric_name,
                        lat,
                        lon,
                    ]
                    idx = list(current_db.created.to_dict().keys())[0]
                    old_created = current_db.created.to_dict()[idx]
                    cn = list(current_db.keys())
                    for i in range(len(new_vals)):
                        names = """UPDATE pycheron SET {cn} = '{new}' where created = '{c}'""".format(
                            cn=cn[i], new=new_vals[i], c=old_created
                        )
                        db.execute(names)

                    self._insert_metric_wrapper(
                        conn, metricStation, timestamp, overwrite=old_created
                    )

                else:
                    print(
                        "Entry already exists, turn on overwrite or change session name"
                    )

                # if there are no gaps, there will be no additional data
                try:
                    for j in range(len(metric[1][i])):
                        metricChannel = metric[1][i][j]

                        snclq = {
                            "network": network,
                            "station": station,
                            "location": location,
                        }

                        metricChannel.update(snclq)

                        timestamp_chan = UTCDateTime(time.time()).isoformat()
                        channel = metricChannel["channel"]
                        metric_name = metricChannel["metric_name"] + "-" + str(j + 1)

                        current_db = self.view(
                            metric_name,
                            network,
                            station,
                            channel,
                            location,
                            session=session,
                            time=(st, et),
                        )

                        # if not currently in db
                        if current_db.empty:
                            names = """INSERT INTO pycheron( network, station, channel, location, session, metric, created,
                             start_time, end_time, latitude, longitude ) VALUES ('{n}','{s}', '{c}', '{l}', '{sn}', '{metric}', '{created}', '{st}', '{et}', '{lat}', '{lon}')""".format(
                                n=network,
                                s=station,
                                c=channel,
                                metric=metric_name,
                                created=timestamp_chan,
                                l=location,
                                sn=self._session_name,
                                st=st,
                                et=et,
                                lat=lat,
                                lon=lon,
                            )
                            db.execute(names)

                            self._insert_metric_wrapper(
                                conn, metricChannel, timestamp_chan
                            )

                        elif self._overwrite:
                            new_vals = [
                                timestamp_chan,
                                network,
                                station,
                                channel,
                                location,
                                session,
                                st,
                                et,
                                metric_name,
                                lat,
                                lon,
                            ]
                            idx = list(current_db.created.to_dict().keys())[0]
                            old_created = current_db.created.to_dict()[idx]
                            cn = list(current_db.keys())

                            for i in range(len(new_vals)):
                                names = """UPDATE pycheron SET {cn} = '{new}' where created = '{c}'""".format(
                                    cn=cn[i], new=new_vals[i], c=old_created
                                )
                                db.execute(names)

                            self._insert_metric_wrapper(
                                conn,
                                metricChannel,
                                timestamp_chan,
                                overwrite=old_created,
                            )

                        else:
                            print(
                                "Entry already exists, turn on overwrite or change session name"
                            )
                except IndexError:
                    continue

        # deals with sohMetric
        elif metric_name == "sohMetric":
            timestamp = UTCDateTime(time.time()).isoformat()
            network = metric["network"]
            station = metric["station"]
            location = metric["location"]
            st = metric["start_time"]
            et = metric["end_time"]
            lat, lon = get_latlon(client, network, station, st, et)

            if location is None or location == "":
                location = "--"

            # check current db
            current_db = self.view(
                metric_name=metric_name,
                network=network,
                station=station,
                channel=None,
                location=location,
                session=session,
                time=(st, et),
            )

            # if not currently in db
            if current_db.empty:
                names = """INSERT INTO pycheron( network, station, channel, location, session ,metric, created, start_time, end_time, latitude, longitude)
                 VALUES ('{n}','{s}', '{chan}', '{l}', '{sn}','{metric}', '{created}',  '{st}', '{et}', '{lat}', '{lon}')""".format(
                    n=network,
                    s=station,
                    chan=None,
                    l=location,
                    sn=self._session_name,
                    metric=metric_name,
                    created=timestamp,
                    st=st,
                    et=et,
                    lat=lat,
                    lon=lon,
                )
                db.execute(names)
                self._insert_metric_wrapper(conn, metric, timestamp)

            elif self._overwrite:
                new_vals = [
                    timestamp,
                    network,
                    station,
                    None,
                    location,
                    session,
                    st,
                    et,
                    metric_name,
                    lat,
                    lon,
                ]
                idx = list(current_db.created.to_dict().keys())[0]
                old_created = current_db.created.to_dict()[idx]
                cn = list(current_db.keys())

                # reverse so created is changed last, since primary key
                cn.reverse()
                new_vals.reverse()
                for i in range(len(new_vals)):
                    names = """UPDATE pycheron SET {cn} = '{new}' where created = '{c}'""".format(
                        cn=cn[i], new=new_vals[i], c=old_created
                    )
                    db.execute(names)

                self._insert_metric_wrapper(
                    conn, metric, timestamp, overwrite=old_created
                )

            else:
                print("Entry already exists, turn on overwrite or change session name")

        # deals with noiseModels (network and station)
        elif "NoiseModel" in metric_name:
            timestamp = UTCDateTime(time.time()).isoformat()
            network = metric["network"]

            if "station" in metric_name:  # stationNoiseModel
                station = metric["station"]
                location = metric["location"]
                lat, lon = get_latlon(client, network, station)
                if location is None or location == "":
                    location = "--"

                current_db = self.view(
                    metric_name=metric_name,
                    network=network,
                    station=station,
                    location=location,
                    session=session,
                )

                # if not currently in db
                if current_db.empty:
                    names = """INSERT INTO pycheron( network, station, location, session ,metric, created, latitude, longitude) VALUES ('{n}','{s}', '{l}', '{sn}',
                         '{metric}', '{created}', '{lat}', '{lon}')""".format(
                        n=network,
                        s=station,
                        l=location,
                        sn=self._session_name,
                        metric=metric_name,
                        created=timestamp,
                        lat=lat,
                        lon=lon,
                    )

                    db.execute(names)

                    self._insert_metric_wrapper(conn, metric, timestamp)

                elif self._overwrite:
                    new_vals = {
                        "created": timestamp,
                        "network": network,
                        "station": station,
                        "location": location,
                        "session": session,
                        "metric": metric_name,
                        "latitude": lat,
                        "longitude": lon,
                    }

                    idx = list(current_db.created.to_dict().keys())[0]
                    old_created = current_db.created.to_dict()[idx]
                    cn = list(new_vals.keys())

                    for j in range(len(new_vals)):
                        names = """UPDATE pycheron SET {cn} = '{new}' where created = '{c}'""".format(
                            cn=cn[j], new=list(new_vals.values())[j], c=old_created
                        )
                        db.execute(names)

                    self._insert_metric_wrapper(
                        conn, metric, timestamp, overwrite=old_created
                    )
                else:
                    print(
                        "Entry already exists, turn on overwrite or change session name"
                    )

                # commit changes
                conn.commit()

            else:  # networkNoiseModel
                current_db = self.view(
                    metric_name=metric_name, network=network, session=self._session_name
                )

                # if not currently in db
                if current_db.empty:

                    names = """INSERT INTO pycheron( network, session , metric, created) VALUES ('{n}', '{sn}',
                         '{metric}', '{created}')""".format(
                        n=network,
                        sn=self._session_name,
                        metric=metric_name,
                        created=timestamp,
                    )

                    db.execute(names)

                    self._insert_metric_wrapper(conn, metric, timestamp)

                elif self._overwrite:

                    new_vals = {
                        "created": timestamp,
                        "network": network,
                        "session": session,
                        "metric": metric_name,
                    }

                    idx = list(current_db.created.to_dict().keys())[0]
                    old_created = current_db.created.to_dict()[idx]
                    cn = list(new_vals.keys())

                    for j in range(len(new_vals)):
                        names = """UPDATE pycheron SET {cn} = '{new}' where created = '{c}'""".format(
                            cn=cn[j], new=list(new_vals.values())[j], c=old_created
                        )

                        db.execute(names)

                    self._insert_metric_wrapper(
                        conn, metric, timestamp, overwrite=old_created
                    )
                else:
                    print(
                        "Entry already exists, turn on overwrite or change session name"
                    )

                # commit changes
                conn.commit()

        # deals with psd plot
        elif metric_name == "psdPlot":
            if isinstance(metric, dict):
                metric = [metric]
            for i in range(len(metric)):
                timestamp = UTCDateTime(time.time()).isoformat()
                metric_name = metric[i]["metric_name"]
                network = metric[i]["network"]
                station = metric[i]["station"]
                channel = metric[i]["channel"]
                location = metric[i]["location"]
                lat, lon = get_latlon(client, network, station)

                if location is None or location == "":
                    location = "--"

                # check current db
                current_db = self.view(
                    metric_name, network, station, channel, location, session=session
                )

                # if not currently in db
                if current_db.empty:
                    names = """INSERT INTO pycheron( network, station, channel, location, session ,metric, created, latitude, longitude) VALUES ('{n}','{s}','{c}', '{l}', '{sn}',
                         '{metric}', '{created}', '{lat}', '{lon}')""".format(
                        n=network,
                        s=station,
                        c=channel,
                        l=location,
                        sn=self._session_name,
                        metric=metric_name,
                        created=timestamp,
                        lat=lat,
                        lon=lon,
                    )

                    db.execute(names)

                    self._insert_metric_wrapper(conn, metric[i], timestamp)

                elif self._overwrite:
                    new_vals = {
                        "created": timestamp,
                        "network": network,
                        "station": station,
                        "channel": channel,
                        "location": location,
                        "session": session,
                        "metric": metric_name,
                        "latitude": lat,
                        "longitude": lon,
                    }

                    idx = list(current_db.created.to_dict().keys())[0]
                    old_created = current_db.created.to_dict()[idx]
                    cn = list(new_vals.keys())

                    for j in range(len(cn)):
                        names = """UPDATE pycheron SET {cn} = '{new}' where created = '{c}'""".format(
                            cn=cn[j], new=list(new_vals.values())[j], c=old_created
                        )
                        db.execute(names)

                    self._insert_metric_wrapper(
                        conn, metric[i], timestamp, overwrite=old_created
                    )
                else:
                    print(
                        "Entry already exists, turn on overwrite or change session name"
                    )

        # deals with everything else
        else:
            for i in range(len(metric)):
                timestamp = UTCDateTime(time.time()).isoformat()
                metric_name = metric[i]["metric_name"]
                snclq = metric[i]["snclq"]
                st = metric[i]["start_time"]
                et = metric[i]["end_time"]

                network, station, channel, location, quality = parse_snclq(snclq)
                lat, lon = get_latlon(client, network, station, st, et)

                if location is None or location == "":
                    location = "--"

                # check current db
                current_db = self.view(
                    metric_name=metric_name,
                    network=network,
                    station=station,
                    channel=channel,
                    location=location,
                    session=session,
                    time=(st, et),
                )

                # if not currently in db
                if current_db.empty:
                    names = """INSERT INTO pycheron( network, station, channel, location, session ,metric, created, start_time, end_time, latitude, longitude) VALUES ('{n}','{s}','{c}', '{l}', '{sn}',
                     '{metric}', '{created}', '{st}', '{et}', '{lat}', '{lon}')""".format(
                        n=network,
                        s=station,
                        c=channel,
                        l=location,
                        sn=self._session_name,
                        metric=metric_name,
                        created=timestamp,
                        st=st,
                        et=et,
                        lat=lat,
                        lon=lon,
                    )
                    db.execute(names)

                    self._insert_metric_wrapper(conn, metric[i], timestamp)

                elif self._overwrite:

                    new_vals = [
                        timestamp,
                        network,
                        station,
                        channel,
                        location,
                        session,
                        st,
                        et,
                        metric_name,
                        lat,
                        lon,
                    ]

                    idx = list(current_db.created.to_dict().keys())[0]
                    old_created = current_db.created.to_dict()[idx]
                    cn = list(current_db.keys())

                    # reverse so created is changed last, since primary key
                    cn.reverse()
                    new_vals.reverse()

                    for j in range(len(new_vals)):
                        # TODO: This is temp until we figure out how to grab the lat/lon and insert it into the db
                        names = """UPDATE pycheron SET {cn} = '{new}' where created = '{c}'""".format(
                            cn=cn[j], new=new_vals[j], c=old_created
                        )
                        db.execute(names)

                    self._insert_metric_wrapper(
                        conn, metric[i], timestamp, overwrite=old_created
                    )
                else:
                    print(
                        "Entry already exists, turn on overwrite or change session name"
                    )

        # commit changes to database
        conn.commit()

    def to_csv(self, out_dir):
        """Converts all tables into CSV files. 1 file per table

        :param out_dir: directory location to save files
        :type out_dir: str

        """
        self._dump_database_to_spreadsheets(self._name, out_dir)

    def to_dataframe(self, table):
        """Converts specific database table into pandas df

        :param table: table name
        :type table: str

        :return: Pandas DataFrame
        :rtype: `pandas.DataFrame`

        """
        return pd.read_sql_table(table)

    def extract_masks(self, mask):
        """Extracts masks. Most masks are dumped into JSON format, this extracxts them and returns them to lists or dicts

        :param mask: mask from database

        :return: extracted mask

        """
        try:
            mask = json.loads(mask)
            for i in range(len(mask)):
                mask[i][0] = json.loads(mask[i][0])
                mask[i][1] = json.loads(mask[i][1])
            return mask

        except IndexError:
            mask = json.loads(mask)
            return mask

    def _extract(self, list):
        """
        Internal function to extract nested list
        :param list: list
        :return: extracted list
        """
        if len(list[0]) == 1:
            new_list = [i[0] for i in list]

        else:
            new_list = [[i] for i in list[0]]
        return new_list

    def _insert_metric_wrapper(self, conn, metric, timestamp, overwrite=None):
        """Internal function wrapper that calls the correct metric insert function

        :param conn: database connection
        :param metric: metric data
        :param timestamp: created timestamp to be used as Primary Key
        :param overwrite: whether to overwrite or not.

        """
        try:
            metric_name = metric["metric_name"]
        except TypeError:
            try:
                metric_name = metric[0]["metric_name"]
            except TypeError:
                metric_name = metric[0][0]["metric_name"]

        if (
            metric_name == "correlationMetric"
            or metric_name == "crossCorrMetric"
            or metric_name == "transferFunctionMetric"
        ):
            self._insert_2_snclq_metric(conn, metric, timestamp, overwrite)

        elif "gapMetric" in metric_name:
            self._insert_gapMetric(conn, metric, timestamp, overwrite)

        elif metric_name == "sohMetric":
            self._insert_sohMetric(conn, metric, timestamp, overwrite)

        elif "NoiseModel" in metric_name:
            self._insert_noiseModel(conn, metric, timestamp, overwrite)

        elif metric_name == "psdPlot":
            self._insert_psdPlots(conn, metric, timestamp, overwrite)

        else:
            self._insert_metric(conn, metric, timestamp, overwrite)
        # commit changes
        conn.commit()

    def _insert_metric(self, conn, metric, timestamp, overwrite=None):
        """Inserts metric into database for non-special case metrics

        :param conn: database connection
        :param metric: metric data
        :param timestamp: timestamp to be used as primary key
        :param overwrite: overwrite term

        """
        db = conn.cursor()

        c = timestamp
        s = metric["snclq"]
        mn = metric["metric_name"]

        if (
            mn == "deadChanMeanMetric"
            or mn == "deadChannelMetric"
            or mn == "deadChanADFMetric"
        ):
            mn = "deadChannel"

        network, station, channel, location, quality = parse_snclq(s)

        if location is None or location == "":
            location = "--"

        if overwrite is None:
            sql = """insert into {mn} ( created, network, station, channel, location, session ) values ( '{c}', '{n}',
             '{s}', '{chan}', '{l}', '{sn}' )""".format(
                mn=mn,
                c=timestamp,
                n=network,
                s=station,
                chan=channel,
                l=location,
                sn=self._session_name,
            )
            db.execute(sql)

            for key, value in metric.items():
                if (
                    key == "masks"
                    or "mask" in key
                    or key == "rep_amp"
                    or key == "spike_times"
                    or key == "dead_channel_times"
                    or key == "dc_offset_times"
                    or key == "detections"
                ):
                    value = json.dumps(value)

                if key == "uncorrected_psds":
                    for i in range(len(value)):
                        value[i][0] = json.dumps(list(value[i][0]))
                        value[i][1] = json.dumps(list(value[i][1]))

                    value = json.dumps(value)

                update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                    mn=mn, key=key.lower(), val=value, c=timestamp
                )
                db.execute(update)

        # overwrite entry
        else:

            for key, value in metric.items():
                # update metric values

                if (
                    key == "masks"
                    or "mask" in key
                    or key == "rep_amp"
                    or key == "spike_times"
                    or key == "dead_channel_times"
                    or key == "dc_offset_times"
                    or key == "detections"
                ):
                    value = json.dumps(value)

                if key == "uncorrected_psds":
                    for i in range(len(value)):
                        value[i][0] = json.dumps(list(value[i][0]))
                        value[i][1] = json.dumps(list(value[i][1]))

                    value = json.dumps(value)

                update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                    mn=mn, key=key, val=value, c=overwrite
                )
                db.execute(update)

            # update session, not found in results dict
            update = (
                """UPDATE {mn} SET session = '{val}' where created = '{c}'""".format(
                    mn=mn, val=self._session_name, c=overwrite
                )
            )
            db.execute(update)

            # update created, must be last
            update = (
                """UPDATE {mn} SET created = '{val}' where created = '{c}'""".format(
                    mn=mn, val=timestamp, c=overwrite
                )
            )
            db.execute(update)

        conn.commit()

    def _insert_2_snclq_metric(self, conn, metric, timestamp, overwrite=None):
        """Internal function which inserts metrics with 2 snclqs (crossCorr, corrMetric, transferFunction)

        :param conn: database conncetion
        :param metric: metric data
        :param timestamp: timestamp to be used as primary key
        :param overwrite: overwrite term
        """
        db = conn.cursor()

        c = timestamp
        s = metric["snclq"]

        snclq1, snclq2 = parse_snclq(s, comb=True)
        network, station, channel, location, quality = parse_snclq(snclq1)
        network2, station2, channel2, location2, quality2 = parse_snclq(snclq2)

        if location == "" and location2 == "":
            location = "--"
            location2 = "--"

        n = "{n}:{n2}".format(n=network, n2=network2)
        s = "{s}:{s2}".format(s=station, s2=station2)
        ch = "{c}:{c2}".format(c=channel, c2=channel2)
        l = "{l}:{l2}".format(l=location, l2=location2)

        mn = metric["metric_name"]

        if overwrite is None:
            sql = """insert into {mn} ( created, network, station, channel, location, session ) values ( '{c}', '{n}',
             '{s}', '{chan}', '{l}', '{sn}' )""".format(
                mn=mn, c=timestamp, n=n, s=s, chan=ch, l=l, sn=self._session_name
            )
            db.execute(sql)

            for key, value in metric.items():

                if key == "masks" or "mask" in key:
                    value = json.dumps(value)

                if mn == "crossCorrMetric" and "starttime" in key:
                    key = "start_time"
                    value = (
                        metric["snclq1_starttime"] + ":" + metric["snclq2_starttime"]
                    )

                elif mn == "crossCorrMetric" and "endtime" in key:
                    key = "end_time"
                    value = metric["snclq1_endtime"] + ":" + metric["snclq2_endtime"]

                update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                    mn=mn, key=key, val=value, c=timestamp
                )
                db.execute(update)

        # overwrite entry
        else:

            for key, value in metric.items():
                # update metric values

                if key == "masks" or "mask" in key:
                    value = json.dumps(value)

                if mn == "crossCorrMetric" and "starttime" in key:
                    key = "start_time"
                    value = (
                        metric["snclq1_starttime"] + ":" + metric["snclq2_starttime"]
                    )

                elif mn == "crossCorrMetric" and "endtime" in key:
                    key = "end_time"
                    value = metric["snclq1_endtime"] + ":" + metric["snclq2_endtime"]

                update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                    mn=mn, key=key, val=value, c=overwrite
                )
                db.execute(update)

            # update session, not found in results dict
            update = (
                """UPDATE {mn} SET session = '{val}' where created = '{c}'""".format(
                    mn=mn, val=self._session_name, c=overwrite
                )
            )
            db.execute(update)

            # update created, must be last
            update = (
                """UPDATE {mn} SET created = '{val}' where created = '{c}'""".format(
                    mn=mn, val=timestamp, c=overwrite
                )
            )
            db.execute(update)

        conn.commit()

    def _insert_gapMetric(self, conn, metric, timestamp, overwrite=None):
        """Internal function which inserts gapMetric into database

        :param conn: database conncetion
        :param metric: gap metric data
        :param timestamp: timestamp to be used as primary key
        :param overwrite: overwrite term
        """
        db = conn.cursor()

        network = metric["network"]
        station = metric["station"]
        location = metric["location"]
        mn = metric["metric_name"]

        if mn == "gapMetric":
            mn = "gapMetricStation"

        if location == None or location == "" or type(location) == str:
            metric["location"] = "--"
            location = "--"

        if overwrite is None:

            sql = """insert into {mn} ( created, network, station, location, session ) values ( '{c}', '{n}',
                     '{s}', '{l}', '{sn}' )""".format(
                mn=mn,
                c=timestamp,
                n=network,
                s=station,
                l=location,
                sn=self._session_name,
            )

            db.execute(sql)

            for key, value in metric.items():

                if key == "masks" or "mask" in key:
                    value = json.dumps(value)

                update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                    mn=mn, key=key, val=value, c=timestamp
                )

                db.execute(update)

        # overwrite entry
        else:

            for key, value in metric.items():
                # update metric values

                if key == "masks" or "mask" in key:
                    value = json.dumps(value)

                update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                    mn=mn, key=key, val=value, c=overwrite
                )
                db.execute(update)

            # update session, not found in results dict
            update = (
                """UPDATE {mn} SET session = '{val}' where created = '{c}'""".format(
                    mn=mn, val=self._session_name, c=overwrite
                )
            )
            db.execute(update)

            # update created, must be last
            update = (
                """UPDATE {mn} SET created = '{val}' where created = '{c}'""".format(
                    mn=mn, val=timestamp, c=overwrite
                )
            )
            db.execute(update)

        conn.commit()

    def _insert_sohMetric(self, conn, metric, timestamp, overwrite=None):
        """Internal function which inserts sohMetric into database

        :param conn: database conncetion
        :param metric: soh metric data
        :param timestamp: timestamp to be used as primary key
        :param overwrite: overwrite term
        """
        db = conn.cursor()

        c = timestamp
        network = metric["network"]
        station = metric["station"]
        location = metric["location"]
        st = metric["start_time"]
        et = metric["end_time"]

        if location == None or location == "":
            location = "--"

        for i in list(metric.keys()):

            if i == "data_quality_flags":
                mn = "sohMetricDataQualityFlags"

                if overwrite is None:
                    sql = """insert into {mn} ( created, network, station, location, session, start_time, end_time ) values ( '{c}', '{n}',
                     '{s}', '{l}', '{sn}', '{st}', '{et}' )""".format(
                        mn=mn,
                        c=timestamp,
                        n=network,
                        s=station,
                        l=location,
                        sn=self._session_name,
                        st=st,
                        et=et,
                    )
                    db.execute(sql)

                    for key, value in metric[i].items():
                        update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                            mn=mn, key=key, val=value, c=timestamp
                        )
                        db.execute(update)

                else:

                    for key, value in metric[i].items():
                        # update metric values

                        update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                            mn=mn, key=key, val=value, c=overwrite
                        )
                        db.execute(update)

                    # update session, not found in results dict
                    update = """UPDATE {mn} SET session = '{val}' where created = '{c}'""".format(
                        mn=mn, val=self._session_name, c=overwrite
                    )
                    db.execute(update)

                    # update created, must be last
                    update = """UPDATE {mn} SET created = '{val}' where created = '{c}'""".format(
                        mn=mn, val=timestamp, c=overwrite
                    )
                    db.execute(update)

            elif i == "activity_flags":

                mn = "sohMetricActivityFlags"

                if overwrite is None:
                    sql = """insert into {mn} ( created, network, station, location, session, start_time, end_time ) values ( '{c}', '{n}',
                     '{s}', '{l}', '{sn}', '{st}', '{et}' )""".format(
                        mn=mn,
                        c=timestamp,
                        n=network,
                        s=station,
                        l=location,
                        sn=self._session_name,
                        st=st,
                        et=et,
                    )
                    db.execute(sql)

                    for key, value in metric[i].items():
                        update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                            mn=mn, key=key, val=value, c=timestamp
                        )
                        db.execute(update)

                else:

                    for key, value in metric[i].items():
                        # update metric values

                        update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                            mn=mn, key=key, val=value, c=overwrite
                        )
                        db.execute(update)

                    # update session, not found in results dict
                    update = """UPDATE {mn} SET session = '{val}' where created = '{c}'""".format(
                        mn=mn, val=self._session_name, c=overwrite
                    )
                    db.execute(update)

                    # update created, must be last
                    update = """UPDATE {mn} SET created = '{val}' where created = '{c}'""".format(
                        mn=mn, val=timestamp, c=overwrite
                    )
                    db.execute(update)

            elif i == "io_clock_flags":

                mn = "sohMetricIOClockFlags"

                if overwrite is None:
                    sql = """insert into {mn} ( created, network, station, location, session, start_time, end_time ) values ( '{c}', '{n}',
                     '{s}', '{l}', '{sn}', '{st}', '{et}' )""".format(
                        mn=mn,
                        c=timestamp,
                        n=network,
                        s=station,
                        l=location,
                        sn=self._session_name,
                        st=st,
                        et=et,
                    )
                    db.execute(sql)

                    for key, value in metric[i].items():
                        update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                            mn=mn, key=key, val=value, c=timestamp
                        )
                        db.execute(update)

                else:

                    for key, value in metric[i].items():
                        # update metric values

                        update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                            mn=mn, key=key, val=value, c=overwrite
                        )
                        db.execute(update)

                    # update session, not found in results dict
                    update = """UPDATE {mn} SET session = '{val}' where created = '{c}'""".format(
                        mn=mn, val=self._session_name, c=overwrite
                    )
                    db.execute(update)

                    # update created, must be last
                    update = """UPDATE {mn} SET created = '{val}' where created = '{c}'""".format(
                        mn=mn, val=timestamp, c=overwrite
                    )
                    db.execute(update)

        conn.commit()

    def _insert_noiseModel(self, conn, metric, timestamp, overwrite=None):
        """Internal function which inserts networkNoiseModel or stationNoiseMNodel into database

        :param conn: database conncetion
        :param metric: metric data
        :param timestamp: timestamp to be used as primary key
        :param overwrite: overwrite term
        """
        db = conn.cursor()

        mn = metric["metric_name"]
        network = metric["network"]

        if overwrite is None:
            if mn == "networkNoiseModel":
                sql = """insert into {mn} ( created, network, session ) values ( '{c}', '{n}', '{sn}' )""".format(
                    mn=mn, c=timestamp, n=network, sn=self._session_name
                )
            else:  # stationNoiseModel
                station = metric["station"]
                location = metric["location"]

                if location == None or location == "":
                    location = "--"

                sql = """insert into {mn} ( created, network, station, location, session ) values ( '{c}', '{n}','{sta}', '{l}','{sn}' )""".format(
                    mn=mn,
                    c=timestamp,
                    n=network,
                    sta=station,
                    l=location,
                    sn=self._session_name,
                )

            db.execute(sql)

            for key, value in metric.items():
                if "stations" in key:
                    value = json.dumps(value)
                if key == "psdStats":
                    for key_, val in value.items():
                        if key_ != "snclq":
                            for j in range(len(value[key_])):
                                for k in range(len(value[key_][j])):
                                    value[key_][j][k] = json.dumps(
                                        list(value[key_][j][k])
                                    )
                        else:
                            for j in range(len(value[key_])):
                                value[key_][j] = json.dumps(list(value[key_][j]))

                    value = json.dumps(value)

                update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                    mn=mn, key=key, val=value, c=timestamp
                )
                db.execute(update)

        else:

            for key, value in metric.items():
                # update metric values

                if "stations" in key or key == "noise_freq_matrix":
                    value = json.dumps(value)

                update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                    mn=mn, key=key, val=value, c=overwrite
                )

                db.execute(update)

            # update session, not found in results dict
            update = (
                """UPDATE {mn} SET session = '{val}' where created = '{c}'""".format(
                    mn=mn, val=self._session_name, c=overwrite
                )
            )
            db.execute(update)

            # update created, must be last
            update = (
                """UPDATE {mn} SET created = '{val}' where created = '{c}'""".format(
                    mn=mn, val=timestamp, c=overwrite
                )
            )
            db.execute(update)

    def _insert_psdPlots(self, conn, metric, timestamp, overwrite=None):
        """Internal function which inserts psdPlots into database

        :param conn: database conncetion
        :param metric: metric data
        :param timestamp: timestamp to be used as primary key
        :param overwrite: overwrite term
        """

        db = conn.cursor()

        mn = metric["metric_name"]
        metric_name = metric["metric_name"]
        network = metric["network"]
        station = metric["station"]
        channel = metric["channel"]
        location = metric["location"]

        if location == None or location == "":
            location = "--"

        if overwrite is None:

            sql = """insert into {mn} ( created, network, station, location, session ) values ( '{c}', '{n}','{sta}', '{l}','{sn}' )""".format(
                mn=mn,
                c=timestamp,
                n=network,
                sta=station,
                l=location,
                sn=self._session_name,
            )

            db.execute(sql)

            for key, value in metric.items():
                if key == "location":
                    if value == None or value == "":
                        value = "--"
                if key == "type":
                    value = json.dumps(value)

                update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                    mn=mn, key=key, val=value, c=timestamp
                )
                db.execute(update)

        else:

            for key, value in metric.items():

                if key == "location":
                    if value == None or value == "":
                        value = "--"

                if key == "type":
                    value = json.dumps(value)

                # update metric values

                update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                    mn=mn, key=key, val=value, c=overwrite
                )

                db.execute(update)

            # update session, not found in results dict
            update = (
                """UPDATE {mn} SET session = '{val}' where created = '{c}'""".format(
                    mn=mn, val=self._session_name, c=overwrite
                )
            )
            db.execute(update)

            # update created, must be last
            update = (
                """UPDATE {mn} SET created = '{val}' where created = '{c}'""".format(
                    mn=mn, val=timestamp, c=overwrite
                )
            )
            db.execute(update)

    # def _insert_calibration(self, conn, metric, timestamp, overwrite=None):
    #     """Internal function which inserts calibration metric into database
    #
    #     :param conn: database conncetion
    #     :param metric: calibration metric metric data
    #     :param timestamp: timestamp to be used as primary key
    #     :param overwrite: overwrite term
    #     """
    #     db = conn.cursor()
    #
    #     for i in range(len(metric[1])):
    #         metric = metric[1][i]
    #         mn = metric[0]["metric_name"]
    #         snclq = metric[0]["snclq"]
    #
    #         network, station, channel, location, quality = parse_snclq(snclq)
    #
    #         if location == None or location == "":
    #             location = "--"
    #
    #         if overwrite is None:
    #             sql = """insert into {mn} ( created, network, station, channel, location, session ) values ( '{c}', '{n}',
    #              '{s}', '{chan}', '{l}', '{sn}' )""".format(mn=mn, c=timestamp, n=network, s=station, chan=channel,
    #                                                         l=location, sn=self._session_name)
    #             db.execute(sql)
    #
    #             for key, value in metric.iteritems():
    #
    #                 if key == "snclq":
    #                     continue
    #
    #                 else:
    #                     update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(mn=mn, key=key,
    #                                                                                                   val=value,
    #                                                                                                   c=timestamp)
    #                     db.execute(update)
    #
    #         # overwrite entry
    #         else:
    #
    #             for key, value in metric.iteritems():
    #                 # update metric values
    #                 if key == "snclq":
    #                     continue
    #                 else:
    #                     update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(mn=mn, key=key,
    #                                                                                                   val=value,
    #                                                                                                   c=overwrite)
    #                     db.execute(update)
    #
    #             # update session, not found in results dict
    #             update = """UPDATE {mn} SET session = '{val}' where created = '{c}'""".format(mn=mn,
    #                                                                                           val=self._session_name,
    #                                                                                           c=overwrite)
    #             db.execute(update)
    #
    #             # update created, must be last
    #             update = """UPDATE {mn} SET created = '{val}' where created = '{c}'""".format(mn=mn, val=timestamp,
    #                                                                                           c=overwrite)
    #             db.execute(update)
    #
    #         conn.commit()
    #
    #     else:
    #         pass

    def _dump_database_to_spreadsheets(self, db_name, filepath):
        """Internal function which converts db into csv files

        :param db_name: database name
        :param filepath: location to save files
        :return:
        """
        db = sqlite3.connect(db_name)
        shortname, extension = os.path.splitext(filepath)
        os.path.isdir(shortname) or os.mkdir(shortname)
        cursor = db.cursor()
        for table in self._list_tables(cursor):
            sheetfile = "%s.csv" % table
            sheetpath = os.path.join(shortname, sheetfile)
            self._dump_table_to_spreadsheet(cursor, table, sheetpath)

    def _list_tables(self, cursor):
        """Internal function used to list tables in database

        :param cursor: database connection
        :return: list of tables
        """
        cursor.execute("select name from sqlite_master")
        return [
            r[0]
            for r in cursor
            if not r[0].startswith("sqlite") and not r[0].endswith("idx")
        ]

    def _dump_table_to_spreadsheet(self, cursor, tablename, sheetpath):
        """Internal function to dump data from tables into spreadsheet

        :param cursor: database connection
        :param tablename: name of table
        :param sheetpath: path to save csv
        :return:
        """
        output = UnicodeWriter(open(sheetpath, "w"))
        cursor.execute("select * from %s" % tablename)
        output.writerow([col[0] for col in cursor.description])
        [_f for _f in (output.writerow(row) for row in cursor) if _f]


class UnicodeWriter:
    """
    A CSV writer which will write rows to CSV file "f",
    which is encoded in the given encoding.
    Source: http://docs.python.org/library/csv.html#csv-examples
    Modified to cope with non-string columns.
    """

    def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
        self.queue = io.StringIO()
        self.writer = csv.writer(self.queue, dialect=dialect, **kwds)
        self.stream = f
        self.encoder = codecs.getincrementalencoder(encoding)()

    def encodeone(self, item):
        if type(item) == str:
            return self.encoder.encode(item)
        else:
            return item

    def writerow(self, row):
        self.writer.writerow([self.encodeone(s) for s in row])
        data = self.queue.getvalue()
        self.stream.write(data)
        self.queue.truncate(0)

    def writerows(self, rows):
        for row in rows:
            self.writerow(row)
