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

from pycheron.util.format import parse_snclq
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
    :param manual: Tell database to use latlon_config.toml for custom lat/lon values
    :type manual: bool
    :param wfdb_conn: Wfdisc database connection
    :type wfdb_conn: WfdiscDB

    To create a database:

    .. code-block:: python

        >>> db = Database("test.db", session_name="session1", overwrite=True)


    """

    def __init__(self, db_name, session_name=None, overwrite=True, manual=False, wfdb_conn=None):
        self._name = db_name
        self._session_name = session_name
        self._overwrite = overwrite
        self._manual = manual
        self._wfdb_conn = wfdb_conn
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
                               min_threshold_exceeded varchar,
                               max_threshold_exceeded varchar,
                               median_threshold_exceeded varchar,
                               mean_threshold_exceeded varchar,
                               variance_threshold_exceeded varchar, 
                               std_threshold_exceeded varchar,
                               highAmp varchar, 
                               pegged varchar,
                               highAmp_rms_mask varchar,
                               pegged_mask varchar,
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

            # make dailyDCOffSetMetricTimes table
            db.execute(
                """create table dailyDCOffset (
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
                               daily_dc_offset_value varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))"""
            )

            # make deadChanADFMetric table
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
                               station_percent_availability_masks varchar,
                               channel_percent_availability_masks varchar,

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
                               dead_channel varchar,
                               low_amp varchar,
                               noise1 varchar,
                               noise2 varchar,
                               highAmp varchar,
                               badResp varchar,
                               dead_channel_gsn_mask varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel))"""
            )

              # make psdMetricInfra table
            db.execute(
                """create table psdMetricInfra(
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
                               percent_above_idc_hnm varchar,
                               dead_channel_linear varchar,
                               uncorrected_psds varchar,
                               percent_below_idc_lnm varchar,
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
                               dead_channel varchar,
                               low_amp varchar,
                               noise1 varchar,
                               noise2 varchar,
                               highAmp varchar,
                               badResp varchar,
                               dead_channel_gsn_mask varchar,
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

                               calibration_signal_counts varchar,
                               event_begin_counts varchar,
                               event_end_counts varchar,
                               event_in_progress_counts varchar,
                               negative_leap_counts varchar,
                               positive_leap_counts varchar,
                               time_correction_applied_counts varchar,

                               calibration_signal_percentages varchar,
                               event_begin_percentages varchar,
                               event_end_percentages varchar,
                               event_in_progress_percentages varchar,
                               negative_leap_percentages varchar,
                               positive_leap_percentages varchar,
                               time_correction_applied_percentages varchar, 

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

                               clock_locked_counts varchar,
                               end_time_series_counts varchar,
                               long_record_read_counts varchar,
                               short_record_read_counts varchar,
                               start_time_series_counts varchar,
                               station_volume_counts varchar,

                               clock_locked_percentages varchar,
                               end_time_series_percentages varchar,
                               long_record_read_percentages varchar,
                               short_record_read_percentages varchar,
                               start_time_series_percentages varchar,
                               station_volume_percentages varchar,

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

                               amplifier_saturation_counts varchar,
                               digital_filter_charging_counts varchar,
                               digitizer_clipping_counts varchar,
                               glitches_counts varchar,
                               missing_padded_data_counts varchar,
                               spikes_counts varchar,
                               suspect_time_tag_counts varchar,
                               telemetry_sync_error_counts varchar,

                               amplifier_saturation_percentages varchar,
                               digital_filter_charging_percentages varchar,
                               digitizer_clipping_percentages varchar,
                               glitches_percentages varchar,
                               missing_padded_data_percentages varchar,
                               spikes_percentages varchar,
                               suspect_time_tag_percentages varchar,
                               telemetry_sync_error_percentages varchar,

                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station))

               """
            )

            # make sohMasksAndGeneral table
            db.execute(
                """create table sohMasksAndGeneral (
                    created VARCHAR PRIMARY KEY,
                    network varchar,
                    station varchar,
                    location varchar,
                    session varchar,
                    start_time varchar,
                    end_time varchar,
                    metric_name varchar,

                    timing_correction varchar,
                    timing_correction_count varchar,
                    timing_quality_record_count varchar,
                    timing_quality_statistics varchar,
                    record_count varchar,
                    num_records_used varchar,
                    noTime varchar,
                    poorTQ varchar,
                    suspectTime varchar,
                    ampSat varchar,
                    digFilterChg varchar,
                    clip varchar,
                    spikes varchar,
                    glitch varchar,
                    missingPad varchar,
                    tsyncErrors varchar,
                    calib varchar,
                    timingCor varchar,
                    noTimeMasks varchar,
                    poorTQMasks varchar,
                    suspectTimeMasks varchar,
                    ampSatMasks varchar,
                    digFilterChgMasks varchar,
                    clipMasks varchar,
                    spikesMasks varchar,
                    glitchMasks varchar,
                    missingPadMasks varchar,
                    tsyncMasks varchar,
                    calibMasks varchar,
                    tcMasks varchar,

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

            # make summaryReportCountsAndValues table
            db.execute(
                """create table summaryReportCountsAndValues(
                                created VARCHAR PRIMARY KEY,
                                range varchar,
                                session varchar,
                                network varchar,
                                station varchar,
                                channel varchar,

                                --basicStatsMetric
                                min_mask int,
                                max_mask int,
                                median_mask int,
                                variance_mask int,
                                std_mask int,
                                rms_mask int,
                                min_threshold_exceeded int,
                                max_threshold_exceeded int,
                                median_threshold_exceeded int,
                                mean_threshold_exceeded int,
                                variance_threshold_exceeded int, 
                                std_threshold_exceeded int,
                                highAmp int, 
                                pegged int,
                                highAmp_rms_mask int,
                                pegged_mask int,

                                --calibrationMetric
                                num_cals_detected int,
                                
                                --dailyDCOffset
                                daily_dc_offset_value int,
                                
                                --DCOffsetTimesMetric
                                dc_offset_times int,

                                --deadChannel
                                dead_chan_mean_metric_masks int,
                                dead_chan_ADF_metric_masks int,
                                is_dead_channel int,

                                --gapMetricStation
                                station_completeness varchar,
                                channel_percent_available varchar,
                                total_gaps int,
                                total_overlaps int,
                                maximum_gap varchar,
                                maximum_overlap varchar,

                                --psdMetric
                                dc_mask int,
                                low_amp_mask int,
                                noise1_mask int,
                                noise2_mask int,
                                hi_amp_mask int,
                                bad_resp_mask int,
                                percent_above_nhnm varchar,
                                percent_below_nlnm varchar,
                                dead_channel_exponent varchar,
                                dead_channel_linear varchar,
                                dead_channel_gsn varchar,
                                --(hourly psd)
                                dead_chan_exp_hourly_masks int,
                                dead_chan_lin_hourly_masks int,
                                dead_chan_gsn_hourly_masks int,

                                --psdMetricInfra
                                dc_mask int,
                                low_amp_mask int,
                                noise1_mask int,
                                noise2_mask int,
                                hi_amp_mask int,
                                bad_resp_mask int,
                                percent_above_idc_hnm varchar,
                                percent_below_idc_lnm varchar,
                                dead_channel_exponent varchar,
                                dead_channel_linear varchar,
                                dead_channel_gsn varchar,
                                --(hourly psd)
                                dead_chan_exp_hourly_masks int,
                                dead_chan_lin_hourly_masks int,
                                dead_chan_gsn_hourly_masks int,

                                --repeatedAmplitudeMetric
                                repAmp int,

                                --spikesMetric
                                total_spike_count int,
                                non_adjacent_spikes varchar,

                                --sohMetricActivityFlags
                                calibration_signal_counts int,
                                event_begin_counts int,
                                event_end_counts int,
                                event_in_progress_counts int,
                                negative_leap_counts int,
                                positive_leap_counts int,
                                time_correction_applied_counts int,

                                calibration_signal_percentages varchar,
                                event_begin_percentages varchar,
                                event_end_percentages varchar,
                                event_in_progress_percentages varchar,
                                negative_leap_percentages varchar,
                                positive_leap_percentages varchar,
                                time_correction_applied_percentages varchar, 
                                
                                --sohMetricDataQualityFlags
                                amplifier_saturation_counts int,
                                digital_filter_charging_counts int,
                                digitizer_clipping_counts int,
                                glitches_counts int,
                                missing_padded_data_counts int,
                                spikes_counts int,
                                suspect_time_tag_counts int,
                                telemetry_sync_error_counts int,

                                amplifier_saturation_percentages varchar,
                                digital_filter_charging_percentages varchar,
                                digitizer_clipping_percentages varchar,
                                glitches_percentages varchar,
                                missing_padded_data_percentages varchar,
                                spikes_percentages varchar,
                                suspect_time_tag_percentages varchar,
                                telemetry_sync_error_percentages varchar,
                                
                                --sohMetricIOClockFlags
                                clock_locked_counts int,
                                end_time_series_counts int,
                                long_record_read_counts int,
                                short_record_read_counts int,
                                start_time_series_counts int,
                                station_volume_counts int,

                                clock_locked_percentages varchar,
                                end_time_series_percentages varchar,
                                long_record_read_percentages varchar,
                                short_record_read_percentages varchar,
                                start_time_series_percentages varchar,
                                station_volume_percentages varchar,

                                --sohMasksAndGeneral
                                timing_correction varchar,
                                timing_correction_count int,
                                timing_quality_record_count varchar,
                                timing_quality_statistics varchar,
                                record_count varchar,
                                num_records_used varchar,
                                noTime varchar,
                                poorTQ varchar,
                                suspectTime varchar,
                                ampSat varchar,
                                digFilterChg varchar,
                                clip varchar,
                                spikes varchar,
                                glitch varchar,
                                missingPad varchar,
                                tsyncErrors varchar,
                                calib varchar,
                                timingCor varchar,
                                noTimeMasks varchar,
                                poorTQMasks varchar,
                                suspectTimeMasks varchar,
                                ampSatMasks varchar,
                                digFilterChgMasks varchar,
                                clipMasks varchar,
                                spikesMasks varchar,
                                glitchMasks varchar,
                                missingPadMasks varchar,
                                tsyncMasks varchar,
                                calibMasks varchar,
                                tcMasks varchar,

                                --dailyPdfPlot
                                noise_masks int,
                                microseism_masks int,
                                banded_masks int,

                                --snrMetric
                                snr_masks int,
                                snr varchar,

                                --correlationMetric
                                correlation_coefficient varchar,
                                p_value varchar,

                                --crossCorrMetric
                                peak_correlation varchar,
                                peak_lag varchar,

                                --transferFunctionMetric
                                gain_ratio varchar,
                                phase_diff varchar,
                                ms_coherence varchar,

                                --staltaMetric
                                max_stalta varchar,
                                event_time varchar,

                                --qcMLMetric
                                dropout_fraction varchar,
                                distinct_values_ratio varchar,
                                packet_time_bandwidth_product varchar,
                                frequency_sigma varchar,
                                discontinuity_max_value varchar,
                                artifacts int,

                                --seedChanSpsCompliance
                                is_chan_sps_Seedcompliant int,
                                does_sps_match_data int,

                                --ChanOrientationCompliance
                                is_chan_orientation_compliant int,

                                --verticalChanOrientationCompliance
                                is_vert_chan_orientation_compliant int,

                                --horzChanOrientationCompliance
                                is_horz_chan_orientation_compliant_tr1 int,
                                is_horz_chan_orientation_compliant_tr2 int,

                                --sampleRateRespVerification
                                sample_rate_resp int,

                                --maxRangeMetric
                                max_range varchar)

                """
            )

            # make transferFunctionMetric table
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

            # make qcMLMetric table
            db.execute(
                """create table qcMLMetric(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               channel varchar,
                               session varchar,
                               metric_name varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               dropout_fraction varchar,
                               distinct_values_ratio varchar,
                               packet_time_bandwidth_product varchar,
                               frequency_sigma varchar,
                               discontinuity_max_value varchar,
                               artifacts varchar,
                               qc_ml_results varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel)
            )"""
            )

            # make basicDBIntegrityMetric table
            db.execute(
                """create table dbIntegrityMetric(
                               created VARCHAR PRIMARY KEY,
                               session varchar,
                               metric_name varchar,
                               integrity_results varchar
            )"""
            )

            # make (for metadataComplianceMetric) seedChanSpsCompliance table
            db.execute(
                """create table seedChanSpsCompliance(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               channel varchar,
                               session varchar,
                               metric_name varchar,
                               metric_subname varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               is_chan_sps_Seedcompliant varchar,
                               does_sps_match_data varchar,
                               FOREIGN KEY (metric_subname) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel)
            )"""
            )

            # make (for metadataComplianceMetric) ChanOrientationCompliance table
            db.execute(
                """create table ChanOrientationCompliance(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               channel varchar,
                               session varchar,
                               metric_name varchar,
                               metric_subname varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               is_chan_orientation_compliant varchar,
                               FOREIGN KEY (metric_subname) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel)
            )"""
            )

            # make (for metadataComplianceMetric) verticalChanOrientationCompliance table
            db.execute(
                """create table verticalChanOrientationCompliance(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               channel varchar,
                               session varchar,
                               metric_name varchar,
                               metric_subname varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               is_vert_chan_orientation_compliant varchar,
                               FOREIGN KEY (metric_subname) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel)
            )"""
            )

            # make (for metadataComplianceMetric) horzChanOrientationCompliance table
            db.execute(
                """create table horzChanOrientationCompliance(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               network_tr2 varchar,
                               station varchar,
                               station_tr2 varchar,
                               location varchar,
                               location_tr2 varchar,
                               channel varchar,
                               channel_tr2 varchar,
                               session varchar,
                               metric_name varchar,
                               metric_subname varchar,
                               snclq varchar,
                               snclq_tr2 varchar,
                               start_time varchar,
                               start_time_tr2 varchar,
                               end_time varchar,
                               end_time_tr2 varchar,
                               is_horz_chan_orientation_compliant_tr1 varchar,
                               is_horz_chan_orientation_compliant_tr2 varchar,
                               FOREIGN KEY (metric_subname) REFERENCES pycheron(metric)
            )"""
            )

            # make (for metadataComplianceMetric) sampleRateRespVerification table
            db.execute(
                """create table sampleRateRespVerification(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               channel varchar,
                               session varchar,
                               metric_name varchar,
                               metric_subname varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               sample_rate_resp varchar,
                               FOREIGN KEY (metric_subname) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel)
            )"""
            )

            # make maxRangeMetric table
            db.execute(
                """create table maxRangeMetric(
                               created VARCHAR PRIMARY KEY,
                               network varchar,
                               station varchar,
                               location varchar,
                               channel varchar,
                               session varchar,
                               metric_name varchar,
                               snclq varchar,
                               start_time varchar,
                               end_time varchar,
                               max_range varchar,
                               FOREIGN KEY (metric_name) REFERENCES pycheron(metric),
                               FOREIGN KEY (network) REFERENCES pycheron(network),
                               FOREIGN KEY (station) REFERENCES pycheron(station),
                               FOREIGN KEY (channel) REFERENCES pycheron(channel)
            )"""
            )

        except (OperationalError, ValueError):
            # if the database already exists, connect
            print("{name} already exists. Connecting to {name}...".format(name=self._name))
            self._conn = self.connect(db_name)

    def set_wfdb_conn(self, wfdb_conn):
        self._wfdb_conn = wfdb_conn

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

    def get_wfdb_conn(self):
        return self._wfdb_conn

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
        """View in the main directory table (pycheron) all of the databases contents.
        Sort by combinations of metric name, network, station, channel, locations, time, and session name.
        If no options are given, it will return everything in the pycheron table.

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
        :param time: Tuple of starttime and endtime strings.
            ex. `time = ("2011-01-01T00:00:00","2011-01-01T23:55:13.190000")`
        :type time: tuple
        :param session: session name
        :type session: str

        :return: Pandas DataFrame
        :rtype: `pandas.DataFrame`
        """
        # if no location, change to --, otherwise it will be unsearchable/filterable
        if location is None or location == "":
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
            tb = pd.read_sql_query("""select {cn} from {tn}""".format(tn="pycheron", cn="*"), self._conn)
        # if gets error, try again
        except DatabaseError:
            tb = pd.read_sql_query("""select {cn} from {tn}""".format(tn="pycheron", cn="*"), self._conn)

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
        print("INSIDE GET METRIC")

        

        
        # gets everything from specific metric table
        tb = pd.read_sql_query("""select {cn} from {tn}""".format(tn=metric_name, cn="*"), self._conn)

        #Try to get location out of dataframe and then use that
        try:
            location = tb["location"][0]
        except KeyError:
            location = None

        # if there is no location, set to --, otherwise it is unsearchable. networkNoiseModel does not have a location
        # so skip
        if (
            (location is None or location == "")
            and metric_name != "networkNoiseModel"
            and metric_name != "stationNoiseModel"
        ):
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

        # loops thru input params
        for key, value in d.items():
            # if value is not set, or key is time, skip
            if value is not None and key != "time":
                # matching input param to database values
                if key == "location" and value != "--":
                    temp = tb.loc[tb[key] == value]
                    if temp.empty:
                        temp = tb.loc[tb[key] == "--"]
                    tb = temp
                else:
                    tb = tb.loc[tb[key] == value]
                print(f"key: {key}, value: {value}, TB: {tb}")
            # if key is time, split into start and end time, then filter.
            elif value is not None and key == "time":
                st = value[0]
                et = value[1]
                tb = tb.loc[tb["start_time"] >= st]
                tb = tb.loc[tb["end_time"] <= et]
                print(f"key: {key}, value: {value}, TB: {tb}")

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
        print(metric_name)
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
                    self._manual,
                    self.get_wfdb_conn(),
                )

                lat2, lon2 = get_latlon(
                    client,
                    network2,
                    station2,
                    metric["snclq2_starttime"],
                    metric["snclq2_endtime"],
                    self._manual,
                    self.get_wfdb_conn(),
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
                    self._manual,
                    self.get_wfdb_conn(),
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
                names = """INSERT INTO pycheron
                ( network, station, channel, location, session, metric, created, 
                  start_time, end_time, latitude, longitude )
                 VALUES ('{n}','{s}', '{c}', '{l}', '{sn}', '{metric}', '{created}', 
                 '{st}', '{et}', '{lat}', '{lon}')""".format(
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

                self._insert_metric_wrapper(conn, metric, timestamp, overwrite=old_created)

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
                lat, lon = get_latlon(
                    client,
                    network,
                    station,
                    st,
                    et,
                    self._manual,
                    self.get_wfdb_conn(),
                )
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
                     start_time, end_time, latitude, longitude ) 
                     VALUES ('{n}','{s}', '{c}', '{l}', '{sn}', '{metric}', '{created}', '{st}', '{et}', 
                     '{lat}', '{lon}')""".format(
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

                    self._insert_metric_wrapper(conn, metricStation, timestamp, overwrite=old_created)

                else:
                    print("Entry already exists, turn on overwrite or change session name")

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
                             start_time, end_time, latitude, longitude ) 
                             VALUES ('{n}','{s}', '{c}', '{l}', '{sn}', '{metric}', 
                             '{created}', '{st}', '{et}', '{lat}', '{lon}')""".format(
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

                            self._insert_metric_wrapper(conn, metricChannel, timestamp_chan)

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
                            print("Entry already exists, turn on overwrite or change session name")
                except IndexError:
                    continue

        # deals with sohMetric
        elif metric_name == "sohMetric":
            for el in metric:
                timestamp = UTCDateTime(time.time()).isoformat()
                network = el["network"]
                station = el["station"]
                location = el["location"]
                st = el["start_time"]
                et = el["end_time"]
                lat, lon = get_latlon(
                    client,
                    network,
                    station,
                    st,
                    et,
                    self._manual,
                    self.get_wfdb_conn(),
                )

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
                    names = """INSERT INTO pycheron( network, station, channel, location, 
                    session , metric, created, start_time, end_time, latitude, longitude)
                    VALUES ('{n}','{s}', '{chan}', '{l}', '{sn}','{metric}', 
                    '{created}',  '{st}', '{et}', '{lat}', '{lon}')""".format(
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
                    self._insert_metric_wrapper(conn, el, timestamp)

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

                    self._insert_metric_wrapper(conn, el, timestamp, overwrite=old_created)

                else:
                    print("Entry already exists, turn on overwrite or change session name")

        # deals with noiseModels (network and station)
        elif "NoiseModel" in metric_name:
            timestamp = UTCDateTime(time.time()).isoformat()
            network = metric["network"]

            if "station" in metric_name:  # stationNoiseModel
                station = metric["station"]
                location = metric["location"]
                lat, lon = get_latlon(
                    client,
                    network,
                    station,
                    manual=self._manual,
                    wfdb_conn=self.get_wfdb_conn(),
                )
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
                    names = """INSERT INTO pycheron( network, station, location, session, 
                        metric, created, latitude, longitude) VALUES ('{n}','{s}', '{l}', '{sn}',
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

                    self._insert_metric_wrapper(conn, metric, timestamp, overwrite=old_created)
                else:
                    print("Entry already exists, turn on overwrite or change session name")

                # commit changes
                conn.commit()

            else:  # networkNoiseModel
                current_db = self.view(metric_name=metric_name, network=network, session=self._session_name)

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

                    self._insert_metric_wrapper(conn, metric, timestamp, overwrite=old_created)
                else:
                    print("Entry already exists, turn on overwrite or change session name")

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
                lat, lon = get_latlon(
                    client,
                    network,
                    station,
                    manual=self._manual,
                    wfdb_conn=self.get_wfdb_conn(),
                )

                if location is None or location == "":
                    location = "--"

                # check current db
                current_db = self.view(metric_name, network, station, channel, location, session=session)

                # if not currently in db
                if current_db.empty:
                    names = """INSERT INTO pycheron( network, station, channel, location, session, 
                            metric, created, latitude, longitude) VALUES ('{n}','{s}','{c}', '{l}','{sn}',
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

                    self._insert_metric_wrapper(conn, metric[i], timestamp, overwrite=old_created)
                else:
                    print("Entry already exists, turn on overwrite or change session name")


        elif metric_name == "dbIntegrityMetric":
            # Because only one db can be checked at a time, we can assume that there will only be result
            timestamp = UTCDateTime(time.time()).isoformat()
            metric_name = metric["metric_name"]

            # check current db
            current_db = self.view(metric_name, None, None, None, None, session=session)

            # if not currently in db
            if current_db.empty:
                names = """INSERT INTO pycheron( network, station, channel, location, session, 
                        metric, created, latitude, longitude) VALUES ('{n}','{s}','{c}', '{l}','{sn}',
                        '{metric}', '{created}', '{lat}', '{lon}')""".format(
                    n=None,
                    s=None,
                    c=None,
                    l=None,
                    sn=self._session_name,
                    metric=metric_name,
                    created=timestamp,
                    lat=None,
                    lon=None,
                )

                db.execute(names)

                self._insert_metric_wrapper(conn, metric, timestamp)

            elif self._overwrite:
                new_vals = {
                    "created": timestamp,
                    "network": None,
                    "station": None,
                    "channel": None,
                    "location": None,
                    "session": session,
                    "metric": metric_name,
                    "latitude": None,
                    "longitude": None,
                }

                idx = list(current_db.created.to_dict().keys())[0]
                old_created = current_db.created.to_dict()[idx]
                cn = list(new_vals.keys())

                for j in range(len(cn)):
                    names = """UPDATE pycheron SET {cn} = '{new}' where created = '{c}'""".format(
                        cn=cn[j], new=list(new_vals.values())[j], c=old_created
                    )
                    db.execute(names)

                self._insert_metric_wrapper(conn, metric[i], timestamp, overwrite=old_created)
            else:
                print("Entry already exists, turn on overwrite or change session name")
        
        # deals with everything else
        else:
            for i in range(len(metric)):
                # TODO T: Ultimately a number of these if statements will need refactoring.
                #         replacing them with helper functions is a potential solution,
                #         (see metadataComplianceMetric as example below)
                timestamp = UTCDateTime(time.time()).isoformat()
                # helper to retrieve correct data for metadataComplianceMetrics
                if metric[i]["metric_name"] == "metadataComplianceMetric":
                    metric_name, network, station, channel, location, st, et, lat, lon = self._metadataComplianceMetric_insert_helper(metric[i], client)

                else:
                    metric_name = metric[i]["metric_name"]
                    snclq = metric[i]["snclq"]
                    st = metric[i]["start_time"]
                    et = metric[i]["end_time"]

                    network, station, channel, location, quality = parse_snclq(snclq)
                    lat, lon = get_latlon(
                        client,
                        network,
                        station,
                        st,
                        et,
                        self._manual,
                        self.get_wfdb_conn(),
                    )

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
                    names = """INSERT INTO pycheron( network, station, channel, location, session, 
                        metric, created, start_time, end_time, latitude, longitude) 
                        VALUES ('{n}','{s}','{c}', '{l}', '{sn}',
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

                    self._insert_metric_wrapper(conn, metric[i], timestamp, overwrite=old_created)
                else:
                    print("Entry already exists, turn on overwrite or change session name")

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

        elif metric_name == "dbIntegrityMetric":
            self._insert_dbIntegrityMetric(conn, metric, timestamp, overwrite)
        
        elif metric_name == "metadataComplianceMetric":
            self._insert_metadataComplianceMetric(conn, metric, timestamp, overwrite)

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

        s = metric["snclq"]
        mn = metric["metric_name"]

        # fields for summary report
        summary_count_dict = {
            "basicStatsMetric": [
                "min_mask",
                "max_mask",
                "median_mask",
                "variance_mask",
                "std_mask",
                "rms_mask",
                "min_threshold_exceeded",
                "max_threshold_exceeded",
                "median_threshold_exceeded",
                "mean_threshold_exceeded",
                "variance_threshold_exceeded",
                "std_threshold_exceeded",
                "highAmp",
                "pegged",
                "highAmp_rms_mask",
                "pegged_mask",
            ],
            "calibrationMetric": ["num_cals_detected"],
            "dailyDCOffset": ["daily_dc_offset_value"],
            "DCOffsetTimesMetric": ["dc_offset_times"],
            "deadChannelMetric": ["is_dead_channel"],
            "deadChanMeanMetric": ["masks"],
            "deadChanADFMetric": ["masks"],
            "psdMetric": [
                "dc_mask",
                "low_amp_mask",
                "noise1_mask",
                "noise2_mask",
                "hi_amp_mask",
                "bad_resp_mask",
                "dead_chan_exp_hourly_masks",
                "dead_chan_lin_hourly_masks",
                "dead_chan_gsn_hourly_masks",
            ],
             "psdMetricInfra": [
                "dc_mask",
                "low_amp_mask",
                "noise1_mask",
                "noise2_mask",
                "hi_amp_mask",
                "bad_resp_mask",
                "dead_chan_exp_hourly_masks",
                "dead_chan_lin_hourly_masks",
                "dead_chan_gsn_hourly_masks",
            ],
            "repeatedAmplitudeMetric": ["count"],
            "spikesMetric": ["total_spike_count"],
            "dailyPdfPlot": ["noise_masks", "microseism_masks", "banded_masks"],
            "snrMetric": ["masks"],
            "qcMLMetric": ["artifacts"],
        }

        summary_value_dict = {
            "psdMetric": [
                "percent_above_nhnm",
                "percent_below_nlnm",
                "dead_channel_exponent",
                "dead_channel_linear",
                "dead_channel_gsn",
            ],
            "psdMetricInfra": [
                "percent_above_idc_hnm",
                "percent_below_idc_lnm",
                "dead_channel_exponent",
                "dead_channel_linear",
                "dead_channel_gsn",
            ],
            "spikesMetric": ["non_adjacent_spikes"],
            "snrMetric": ["SNR"],
            "staltaMetric": ["max_stalta", "event_time"],
            "qcMLMetric": [
                "dropout_fraction",
                "distinct_values_ratio",
                "packet_time_bandwidth_product",
                "frequency_sigma",
                "discontinuity_max_value"
            ],
            "maxRangeMetric": ["max_range"],
        }

        if mn == "deadChanMeanMetric" or mn == "deadChannelMetric" or mn == "deadChanADFMetric":
            mn = "deadChannel"

        network, station, channel, location, quality = parse_snclq(s)

        dt_range = None
        if "start_time" and "end_time" in metric:
            dt_range = uniform_time(metric["start_time"], metric["end_time"])
            qc_insert = summary_table_check_exist(db, self._session_name, network, station, channel, dt_range)

            if qc_insert:
                db.execute(qc_insert)

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
                valuej = None
                if (
                    key == "masks"
                    or "mask" in key
                    or key == "rep_amp"
                    or key == "spike_times"
                    or key == "dead_channel_times"
                    or key == "dc_offset_times"
                    or key == "detections"
                    or key == "daily_dc_offset_value"
                    or key == "qc_ml_results"
                ):
                    valuej = json.dumps(value)

                if key == "daily_dc_offset_value" and value:
                    value = [value]

                if key == "uncorrected_psds":
                    for i in range(len(value)):
                        value[i][0] = json.dumps(list(value[i][0]))
                        value[i][1] = json.dumps(list(value[i][1]))

                    valuej = json.dumps(value)

                if valuej:
                    update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                        mn=mn, key=key.lower(), val=valuej, c=timestamp
                    )
                else:
                    update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                        mn=mn, key=key.lower(), val=value, c=timestamp
                    )
                db.execute(update)

                sum_counts_sql = self._insert_summary_sql(mn, key, value, summary_count_dict, dt_range, network, station, "_counts", channel)
                if sum_counts_sql:
                    db.execute(sum_counts_sql)

                sum_value_sql = self._insert_summary_sql(mn, key, value, summary_value_dict, dt_range, network, station, "_value", channel)
                if sum_value_sql:
                    db.execute(sum_value_sql)

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
                    or key == "daily_dc_offset_value"
                    or key == "qc_ml_results"
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

    def _insert_2_snclq_metric(self, conn, metric, timestamp, overwrite=None):
        """Internal function which inserts metrics with 2 snclqs (crossCorr, corrMetric, transferFunction)

        :param conn: database conncetion
        :param metric: metric data
        :param timestamp: timestamp to be used as primary key
        :param overwrite: overwrite term
        """
        db = conn.cursor()
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

        summary_value_dict = {
            "correlationMetric": ["correlation_coefficient", "p_value"],
            "crossCorrMetric": ["peak_correlation", "peak_lag"],
            "transferFunctionMetric": ["gain_ratio", "phase_diff", "ms_coherence"],
        }

        dt_range = None
        if "start_time" and "end_time" in metric:
            dt_range = uniform_time(metric["start_time"], metric["end_time"])
            qc_insert = summary_table_check_exist(db, self._session_name, network, station, channel, dt_range)
            if qc_insert:
                db.execute(qc_insert)
        if "snclq1_starttime" and "snclq1_endtime" in metric:
            dt_range = uniform_time(metric["snclq1_starttime"], metric["snclq1_endtime"])
            # TODO: Add logic for snclq2 if crossCorr for insert query
            qc_insert = summary_table_check_exist(db, self._session_name, network, station, channel, dt_range)
            if qc_insert:
                db.execute(qc_insert)

        if overwrite is None:
            sql = """insert into {mn} ( created, network, station, channel, location, session ) values ( '{c}', '{n}',
             '{s}', '{chan}', '{l}', '{sn}' )""".format(
                mn=mn, c=timestamp, n=n, s=s, chan=ch, l=l, sn=self._session_name
            )
            db.execute(sql)

            for key, value in metric.items():
                valuej = None
                if key == "masks" or "mask" in key:
                    valuej = json.dumps(value)

                if mn == "crossCorrMetric" and "starttime" in key:
                    key = "start_time"
                    value = metric["snclq1_starttime"] + ":" + metric["snclq2_starttime"]

                elif mn == "crossCorrMetric" and "endtime" in key:
                    key = "end_time"
                    value = metric["snclq1_endtime"] + ":" + metric["snclq2_endtime"]

                if valuej:
                    update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                        mn=mn, key=key, val=valuej, c=timestamp
                    )
                else:
                    update = """UPDATE {mn} SET '{key}' = '{val}' where created = '{c}'""".format(
                        mn=mn, key=key, val=value, c=timestamp
                    )
                db.execute(update)

                sum_value_sql = self._insert_summary_sql(mn, key, value, summary_value_dict, dt_range, network, station, "_value", channel)
                if sum_value_sql:
                    db.execute(sum_value_sql)

        # overwrite entry
        else:

            for key, value in metric.items():
                # update metric values

                if key == "masks" or "mask" in key:
                    value = json.dumps(value)

                if mn == "crossCorrMetric" and "starttime" in key:
                    key = "start_time"
                    value = metric["snclq1_starttime"] + ":" + metric["snclq2_starttime"]

                elif mn == "crossCorrMetric" and "endtime" in key:
                    key = "end_time"
                    value = metric["snclq1_endtime"] + ":" + metric["snclq2_endtime"]

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

        summary_count_dict = {
            "gapMetricStation": ["total_gaps", "total_overlaps"],
        }

        summary_value_dict = {
            "gapMetricStation": [
                "station_completeness",
                "channel_percent_available",
                "maximum_gap",
                "maximum_overlap",
            ],
        }

        channel = None
        dt_range = None
        if "start_time" and "end_time" in metric:
            dt_range = uniform_time(metric["start_time"], metric["end_time"])
            qc_insert = summary_table_check_exist(db, self._session_name, network, station, channel, dt_range)

            if qc_insert:
                db.execute(qc_insert)

        if mn == "gapMetric":
            mn = "gapMetricStation"

        if location is None or location == "" or type(location) == str:
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

                sum_counts_sql = self._insert_summary_sql(mn, key, value, summary_count_dict, dt_range, network, station, "_counts", channel)
                if sum_counts_sql:
                    db.execute(sum_counts_sql)

                sum_value_sql = self._insert_summary_sql(mn, key, value, summary_value_dict, dt_range, network, station, "_value", channel)
                if sum_value_sql:
                    db.execute(sum_value_sql)

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

    def _insert_sohMetric(self, conn, metric, timestamp, overwrite=None):
        """Internal function which inserts sohMetric into database

        :param conn: database conncetion
        :param metric: soh metric data
        :param timestamp: timestamp to be used as primary key
        :param overwrite: overwrite term
        """
        db = conn.cursor()

        network = metric["network"]
        station = metric["station"]
        location = metric["location"]
        st = metric["start_time"]
        et = metric["end_time"]

        summary_count_dict = {
            "sohMetric": [
                "calibration_signal_counts",
                "event_begin_counts",
                "event_end_counts",
                "event_in_progress_counts",
                "negative_leap_counts",
                "positive_leap_counts",
                "time_correction_applied_counts",
                "amplifier_saturation_counts",
                "digital_filter_charging_counts",
                "digitizer_clipping_counts",
                "glitches_counts",
                "missing_padded_data_counts",
                "spikes_counts",
                "suspect_time_tag_counts",
                "telemetry_sync_error_counts",
                "clock_locked_counts",
                "end_time_series_counts",
                "long_record_read_counts",
                "short_record_read_counts",
                "start_time_series_counts",
                "station_volume_counts",
                "timing_correction_count",
            ]
        }

        summary_value_dict = {
            "sohMetric": [
                "calibration_signal_percentages",
                "event_begin_percentages",
                "event_end_percentages",
                "event_in_progress_percentages",
                "negative_leap_percentages",
                "positive_leap_percentages",
                "time_correction_applied_percentages",
                "amplifier_saturation_percentages",
                "digital_filter_charging_percentages",
                "digitizer_clipping_percentages",
                "glitches_percentages",
                "missing_padded_data_percentages",
                "spikes_percentages",
                "suspect_time_tag_percentages",
                "telemetry_sync_error_percentages",
                "clock_locked_percentages",
                "end_time_series_percentages",
                "long_record_read_percentages",
                "short_record_read_percentages",
                "start_time_series_percentages",
                "station_volume_percentages",
                "timing_correction",
                "timing_quality_record_count",
                "timing_quality_statistics",
                "record_count",
                "num_records_used",
                "noTime",
                "poorTQ",
                "suspectTime",
                "ampSat",
                "digFilterChg",
                "clip",
                "spikes",
                "glitch",
                "missingPad",
                "tsyncErrors",
                "calib",
                "timingCor",
                "noTimeMasks",
                "poorTQMasks",
                "suspectTimeMasks",
                "ampSatMasks",
                "digFilterChgMasks",
                "clipMasks",
                "spikesMasks",
                "glitchMasks",
                "missingPadMasks",
                "tsyncMasks",
                "calibMasks",
                "tcMasks",
            ],
        }

        channel = None

        dt_range = None
        if "start_time" and "end_time" in metric:
            dt_range = uniform_time(metric["start_time"], metric["end_time"])
            qc_insert = summary_table_check_exist(db, self._session_name, network, station, channel, dt_range)

            if qc_insert:
                db.execute(qc_insert)

        if location is None or location == "":
            location = "--"

        soh_branches = {
            "data_quality": {
                "counts": "data_quality_flag_counts",
                "percentages": "data_quality_percentages",
                "table": "sohMetricDataQualityFlags",
            },
            "io_clock": {
                "counts": "io_clock_flag_counts",
                "percentages": "io_clock_flags_percentages",
                "table": "sohMetricIOClockFlags",
            },
            "activity": {
                "counts": "activity_flag_counts",
                "percentages": "activity_flag_percentages",
                "table": "sohMetricActivityFlags",
            },
        }
        has_inserted = False

        for flag in list(metric.keys()):
            root_key = pull_root_key(flag)

            if root_key:
                count_dict = metric[soh_branches[root_key]["counts"]]
                percent_dict = metric[soh_branches[root_key]["percentages"]]
                mn = soh_branches[root_key]["table"]

                count_key_lst = list(count_dict.keys())
                percent_key_lst = list(percent_dict.keys())

                count_key_appended = [x + "_counts" for x in count_key_lst]
                percent_key_appended = [x + "_percentages" for x in percent_key_lst]

                count_key_str = str(count_key_appended)[1:-1].replace("'", "")
                percent_key_str = str(percent_key_appended)[1:-1].replace("'", "")

                count_values_lst = list(count_dict.values())
                percent_values_lst = list(percent_dict.values())

                count_values_str = str(list(map(lambda x: str(x), count_values_lst)))[1:-1]
                percent_values_str = str(list(map(lambda x: str(x), percent_values_lst)))[1:-1]

                timestamp = UTCDateTime(time.time()).isoformat()

                if overwrite is None:
                    sql = """insert into {mn} ( created, network, station, location, session, 
                             start_time, end_time, metric_name, {count_key_str}, {percent_key_str} ) 
                             values ( '{c}', '{n}','{s}', '{l}', '{sn}', '{st}', '{et}', 
                             '{metric_name}', {count_values_str}, {percent_values_str})""".format(
                        mn=mn,
                        count_key_str=count_key_str,
                        percent_key_str=percent_key_str,
                        c=timestamp,
                        n=network,
                        s=station,
                        l=location,
                        sn=self._session_name,
                        st=st,
                        et=et,
                        metric_name=metric["metric_name"],
                        count_values_str=count_values_str,
                        percent_values_str=percent_values_str,
                    )
                    db.execute(sql)

                    # grab updated metric names, zip them together into a single dict representing percents and counts
                    app = dict(zip(count_key_appended, count_values_lst))
                    app.update(zip(percent_key_appended, percent_values_lst))

                    # grab all other metrics make them a flattened dict with counts and percents to iterate through
                    metric_summary_prep = {k: v for k, v in metric.items() if not isinstance(v, dict)}
                    metric_summary_prep.update(app)
                    metric_summary_prep.update(metric["timing_quality"])
                    mspn = metric_summary_prep["metric_name"]

                    for key, value in metric_summary_prep.items():

                        sum_counts_sql = self._insert_summary_sql(mspn, key, value, summary_count_dict, dt_range, network, station, "_counts", channel)
                        if sum_counts_sql:
                            db.execute(sum_counts_sql)

                        sum_value_sql = self._insert_summary_sql(mspn, key, value, summary_value_dict, dt_range, network, station, "_value", channel)
                        if sum_value_sql:
                            db.execute(sum_value_sql)



                else:
                    for count_key in count_key_lst:
                        count_sql = (
                            """UPDATE {mn} SET '{count_key}' = '{count_values_str}' WHERE created = '{c}'""".format(
                                mn=mn,
                                count_key=count_key + "_counts",
                                count_values_str=count_dict[count_key],
                                c=overwrite,
                            )
                        )
                        db.execute(count_sql)
                    for per_key in percent_key_lst:
                        percent_sql = (
                            """UPDATE {mn} SET '{percent_key}' = '{percent_values_str}' WHERE created = '{c}'""".format(
                                mn=mn,
                                percent_key=per_key + "_counts",
                                percent_values_str=percent_dict[per_key],
                                c=overwrite,
                            )
                        )
                        db.execute(percent_sql)

            # elif here is to account for metrics which aren't counts or percentages, but should still be recorded
            elif (
                flag
                not in (
                    "metric_name",
                    "network",
                    "station",
                    "channel",
                    "location",
                    "start_time",
                    "end_time",
                )
                and not flag.endswith("percentages")
            ):
                mn = "sohMasksAndGeneral"

                if isinstance(metric[flag], dict):
                    for key, value in metric[flag].items():
                        # update metric values
                        if overwrite is None:
                            if has_inserted:
                                update = """Update {mn} SET '{key}' = '{val}' WHERE created = '{c}'""".format(
                                    mn=mn, c=timestamp, key=key, val=value
                                )
                            else:
                                update = """insert into {mn} (created, network, station, location, session, 
                                            start_time, end_time, metric_name, {key}) values ('{c}', '{n}',
                                            '{s}', '{l}', '{sn}', '{st}', '{et}', '{metric_name}', 
                                            '{val}')""".format(
                                    mn=mn,
                                    c=timestamp,
                                    n=network,
                                    s=station,
                                    l=location,
                                    sn=self._session_name,
                                    st=st,
                                    et=et,
                                    metric_name=metric["metric_name"],
                                    key=key,
                                    val=value,
                                )
                                has_inserted = True
                        else:
                            update = """Update {mn} SET '{key}' = '{val}' WHERE created = '{c}'""".format(
                                mn=mn, c=overwrite, key=key, val=value
                            )
                        db.execute(update)

                else:
                    if overwrite is None:
                        if has_inserted:
                            update = """Update {mn} SET '{key}' = '{val}' WHERE created = '{c}'""".format(
                                mn=mn, c=timestamp, key=flag, val=metric[flag]
                            )
                        else:
                            update = """insert into {mn} (created, network, station, location, session, 
                                        start_time, end_time, metric_name, {key}) values ('{c}', '{n}',
                                        '{s}', '{l}', '{sn}', '{st}', '{et}', '{metric_name}', 
                                        '{val}')""".format(
                                mn=mn,
                                c=timestamp,
                                n=network,
                                s=station,
                                l=location,
                                sn=self._session_name,
                                st=st,
                                et=et,
                                metric_name=metric["metric_name"],
                                key=key,
                                val=value,
                            )
                            has_inserted = True

                    else:
                        update = """Update {mn} SET '{key}' = '{val}' WHERE created = '{c}'""".format(
                            mn=mn, c=overwrite, key=flag, val=metric[flag]
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

                if location is None or location == "":
                    location = "--"

                sql = """insert into {mn} ( created, network, station, location, session )
                         values ( '{c}', '{n}','{sta}', '{l}','{sn}' )""".format(
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
                                    value[key_][j][k] = json.dumps(list(value[key_][j][k]))
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
            update = """UPDATE {mn} SET session = '{val}' where created = '{c}'""".format(
                mn=mn, val=self._session_name, c=overwrite
            )
            db.execute(update)

            # update created, must be last
            update = """UPDATE {mn} SET created = '{val}' where created = '{c}'""".format(
                mn=mn, val=timestamp, c=overwrite
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
        network = metric["network"]
        station = metric["station"]
        location = metric["location"]

        if location is None or location == "":
            location = "--"

        if overwrite is None:

            sql = """insert into {mn} ( created, network, station, location, session )
                     values ( '{c}', '{n}','{sta}', '{l}','{sn}' )""".format(
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
                    if value is None or value == "":
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
                    if value is None or value == "":
                        value = "--"

                if key == "type":
                    value = json.dumps(value)

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

    def _insert_dbIntegrityMetric(self, conn, metric, timestamp, overwrite=None):
        db = conn.cursor()
        mn = metric["metric_name"]
        ir = json.dumps(metric["integrity_results"])

        if overwrite is None:
            sql = """insert into {mn} (created, session, metric_name, integrity_results) values ('{c}', '{sn}',
            '{mn}', '{ir}')""".format(
                mn=mn,
                c=timestamp,
                sn=self._session_name,
                ir=ir
            )
            db.execute(sql)
        else:
            sql = """UPDATE {mn} SET integrity_results = '{ir}'""".format(
                mn=mn,
                ir=ir
            )
            db.execute(sql)
        conn.commit()

    def _insert_metadataComplianceMetric(self, conn, metric, timestamp, overwrite=None):
        db = conn.cursor()
        mn = metric["metric_name"]
        smn = metric["metric_subname"]

        network, station, channel, location, _ = parse_snclq(metric["snclq"])

        if location == "":
            location = "--"

            

        # fields for summary report metado
        summary_count_dict = {
            "seedChanSpsCompliance": [
                "is_chan_sps_Seedcompliant",
                "does_sps_match_data",
            ],
            "ChanOrientationCompliance": ["is_chan_orientation_compliant"],
            "verticalChanOrientationCompliance": ["is_vert_chan_orientation_compliant"],
            "horzChanOrientationCompliance": [
                "is_horz_chan_orientation_compliant_tr1",
                "is_horz_chan_orientation_compliant_tr2",
            ],
            "sampleRateRespVerification": ["sample_rate_resp"],
        }
        
        if smn == "seedChanSpsCompliance":
            sql = self._insert_seedChanSpsCompliance_sql(metric, network, station, channel, location, timestamp)
        elif smn == "ChanOrientationCompliance":
            sql = self._insert_ChanOrientationCompliance_sql(metric, network, station, channel, location, timestamp)
        elif smn == "verticalChanOrientationCompliance":
            sql = self._insert_verticalChanOrientationCompliance_sql(metric, network, station, channel, location, timestamp)
        elif smn == "horzChanOrientationCompliance":
            sql = self._insert_horzChanOrientationCompliance_sql(metric, timestamp)
        elif smn == "sampleRateRespVerification":
            sql = self._insert_sampleRateRespVerification_sql(metric, network, station, channel, location, timestamp)
        else:
            raise ValueError(f"{smn} is not a valid metric_subname for {mn}, unable to write to sqlite")

        dt_range=None
        if "start_time" and "end_time" in metric:
            dt_range = uniform_time(metric["start_time"], metric["end_time"])
            qc_insert = summary_table_check_exist(db, self._session_name, network, station, channel, dt_range)

            if qc_insert:
                db.execute(qc_insert)

        db.execute(sql)

        for key, value in metric.items():
            # TODO T: Need case for horzChan and other metrics which have two snclq values
            sum_sql = self._insert_summary_sql(smn, key, value, summary_count_dict, dt_range, network, station, "_counts", channel)
            if sum_sql:
                db.execute(sum_sql)

        conn.commit()

    def _insert_summary_sql(self, metric_name, key, value, summary_dict, dt_range, network, station, c_or_v, channel=None):
        # Determine if counts are values are to be inserted
        if c_or_v not in ("_counts", "_value"):
            raise ValueError(f"_insert_summary_sql function error: Can not insert to summaryCountsAndValues, expect '_counts' or '_value' for c_or_v, but received: {c_or_v}")
        else:
            # If counts, use the counting function to assign them. If values, use the standalone function
            if c_or_v == "_counts":
                assign_to_report = count_included
            if c_or_v == "_value":
                assign_to_report = standalone
        # Iterate through metric keys. If a key is in the summary dictionary, add it to the summary report with appropriate function
        if metric_name in summary_dict.keys():
            if key in summary_dict[metric_name]:
                key = table_key_qc(key, metric_name)
                qc_layout = assign_to_report(key, value)
                update = """UPDATE summaryReportCountsAndValues SET '{qc_layout_key}' = '{qc_layout_value}' 
                    WHERE range = '{dt}' AND session = '{ses}' AND network = '{net}' 
                    AND station = '{sta}'""".format(
                    dt=dt_range,
                    qc_layout_key=key,
                    qc_layout_value=qc_layout[f"{key}{c_or_v}"],
                    ses=self._session_name,
                    net=network,
                    sta=station,
                )
                if channel:
                    update += """AND channel = '{chan}'""".format(chan=channel)
                return update

    def _insert_seedChanSpsCompliance_sql(self, metric, network, station, channel, location, timestamp):
        sql = """insert into {ms} (created, network, station, location, channel, session, metric_name, metric_subname, snclq, start_time, end_time, is_chan_sps_Seedcompliant, does_sps_match_data)
        values ('{c}', '{n}', '{s}', '{l}', '{ch}', '{sn}', '{mn}', '{ms}', '{snclq}', '{st}', '{et}', '{is_chan_sps}', '{does_sps_match}')""".format(
            ms=metric["metric_subname"],
            c=timestamp,
            n=network,
            s=station,
            l=location,
            ch=channel,
            sn=self._session_name,
            mn=metric["metric_name"],
            snclq=metric["snclq"],
            st=metric["start_time"],
            et=metric["end_time"],
            is_chan_sps=metric["is_chan_sps_Seedcompliant"],
            does_sps_match=metric["does_sps_match_data"],
        )
        return sql

    def _insert_ChanOrientationCompliance_sql(self, metric, network, station, channel, location, timestamp):
        sql = """insert into {ms} (created, network, station, location, channel, session, metric_name, metric_subname, snclq, start_time, end_time, is_chan_orientation_compliant)
        values ('{c}', '{n}', '{s}', '{l}', '{ch}', '{sn}', '{mn}', '{ms}', '{snclq}', '{st}', '{et}', '{is_chan_ori}')""".format(
            ms=metric["metric_subname"],
            c=timestamp,
            n=network,
            s=station,
            l=location,
            ch=channel,
            sn=self._session_name,
            mn=metric["metric_name"],
            snclq=metric["snclq"],
            st=metric["start_time"],
            et=metric["end_time"],
            is_chan_ori=metric["is_chan_orientation_compliant"],
        )
        return sql

    def _insert_verticalChanOrientationCompliance_sql(self, metric, network, station, channel, location, timestamp):
        sql = """insert into {ms} (created, network, station, location, channel, session, metric_name, metric_subname, snclq, start_time, end_time, is_vert_chan_orientation_compliant)
        values ('{c}', '{n}', '{s}', '{l}', '{ch}', '{sn}', '{mn}', '{ms}', '{snclq}', '{st}', '{et}', '{is_vchan_ori}')""".format(
            ms=metric["metric_subname"],
            c=timestamp,
            n=network,
            s=station,
            l=location,
            ch=channel,
            sn=self._session_name,
            mn=metric["metric_name"],
            snclq=metric["snclq"],
            st=metric["start_time"],
            et=metric["end_time"],
            is_vchan_ori=metric["is_vert_chan_orientation_compliant"],
        )
        return sql

    def _insert_horzChanOrientationCompliance_sql(self, metric, timestamp):
        # TODO T: Update to be similar to other multi network/station metrics
        # (crossCor, transferFunction, etc.)
        network1, station1, channel1, location1, _ = parse_snclq(metric["snclq"])
        network2, station2, channel2, location2, _ = parse_snclq(metric["snclq_tr2"])

        if location1 == "":
            location1 = "--"
        if location2 == "":
            location2 = "--"

        sql = """insert into {ms} (created, network, network_tr2, station, station_tr2, location,
        location_tr2, channel, channel_tr2, session, metric_name, metric_subname, snclq, snclq_tr2,
        start_time, start_time_tr2, end_time, end_time_tr2, is_horz_chan_orientation_compliant_tr1,
        is_horz_chan_orientation_compliant_tr2) values ('{c}', '{n1}', '{n2}', '{s1}', '{s2}', '{l1}', '{l2}',
        '{ch1}', '{ch2}', '{sn}', '{mn}', '{ms}', '{snclq1}', '{snclq2}', '{st1}', '{st2}', '{et1}', '{et2}',
        '{is_hchan_ori1}', '{is_hchan_ori2}')""".format(
            ms=metric["metric_subname"],
            c=timestamp,
            n1=network1,
            n2=network2,
            s1=station1,
            s2=station2,
            l1=location1,
            l2=location2,
            ch1=channel1,
            ch2=channel2,
            sn=self._session_name,
            mn=metric["metric_name"],
            snclq1=metric["snclq"],
            snclq2=metric["snclq_tr2"],
            st1=metric["start_time"],
            st2=metric["start_time_tr2"],
            et1=metric["end_time"],
            et2=metric["end_time_tr2"],
            is_hchan_ori1=metric["is_horz_chan_orientation_compliant_tr1"],
            is_hchan_ori2=metric["is_horz_chan_orientation_compliant_tr2"],
        )
        return sql

    def _insert_sampleRateRespVerification_sql(self, metric, network, station, channel, location, timestamp):
        sql = """insert into {ms} (created, network, station, location, channel, session, metric_name, metric_subname, snclq, start_time, end_time, sample_rate_resp)
        values ('{c}', '{n}', '{s}', '{l}', '{ch}', '{sn}', '{mn}', '{ms}', '{snclq}', '{st}', '{et}', '{samp_rate_resp}')""".format(
            ms=metric["metric_subname"],
            c=timestamp,
            n=network,
            s=station,
            l=location,
            ch=channel,
            sn=self._session_name,
            mn=metric["metric_name"],
            snclq=metric["snclq"],
            st=metric["start_time"],
            et=metric["end_time"],
            samp_rate_resp=metric["sample_rate_resp"],
        )
        return sql

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
        return [r[0] for r in cursor if not r[0].startswith("sqlite") and not r[0].endswith("idx")]

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

    def _metadataComplianceMetric_insert_helper(self, metric, client):
        metric_name = metric["metric_subname"]

        if metric_name == "horzChanOrientationCompliance":
            network1, station1, channel1, location1, _ = parse_snclq(metric["snclq"])
            network2, station2, channel2, location2, _ = parse_snclq(metric["snclq_tr2"])

            st1 = metric["start_time"]
            st2 = metric["start_time_tr2"]
            et1 = metric["end_time"]
            et2 = metric["end_time_tr2"]

            lat1, lon1 = get_latlon(
                client,
                network1,
                station1,
                st1,
                et1,
                self._manual,
                self.get_wfdb_conn(),
            )

            lat2, lon2 = get_latlon(
                client,
                network2,
                station2,
                st2,
                et2,
                self._manual,
                self.get_wfdb_conn(),
            )

            if location1 is None or location1 == "":
                location1 = "--"
            if location2 is None or location2 == "":
                location2 = "--"
            
            st = f"{metric['start_time']}:{metric['start_time_tr2']}"
            et = f"{metric['end_time']}:{metric['end_time_tr2']}"

            network = f"{network1}:{network2}"
            station = f"{station1}:{station2}"
            channel = f"{channel1}:{channel2}"
            location = f"{location1}:{location2}"
            lat = f"{lat1}:{lat2}"
            lon = f"{lon1}:{lon2}"

            return (metric_name, network, station, channel, location, st, et, lat, lon)
        
        else:
            network, station, channel, location, _ = parse_snclq(metric["snclq"])

            st = metric["start_time"]
            et = metric["end_time"]
            
            lat, lon = get_latlon(
                client,
                network,
                station,
                st,
                et,
                self._manual,
                self.get_wfdb_conn(),
            )

            if location is None or location == "":
                location = "--"

            return (metric_name, network, station, channel, location, st, et, lat, lon)


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


def count_included(key, value):
    """
    Reads in key/value pairs for metrics requiring counts in summaryReportCountsAndValues

    :param key: metric value title
    :type key: str
    :param value: metric value
    :type value: anything but str (excluding string representations of bools/Nones)
    """
    qc_layout = {f"{key}_counts": 0}
    offended = True
    dti = is_datetime(value)
    try:
        if offended:
            qc_layout, offended = count_number(qc_layout, key, value, offended)
        if offended:
            qc_layout, offended = count_bool(qc_layout, key, value, offended)
        if offended:
            qc_layout, offended = count_none(qc_layout, key, value, offended)
        if offended:
            qc_layout, offended = count_datetime(qc_layout, key, dti, offended)
        if offended:
            qc_layout, offended = count_list(qc_layout, key, value, offended)
        if offended:
            qc_layout, offended = count_dict(qc_layout, key, value, offended)
        if offended:
            raise ValueError(
                f"Offending value: {value} --  Type {type(value)} not considered in count_included. -- Key {key}"
            )
    except Exception as e:
        raise e
    return qc_layout


def count_bool(qc_layout, key, value, offended):
    """
    Counts boolean values. If True, 1 count is added. If false, count is 0.

    :param qc_layout: Dict containing count of a given metric
    :type qc_layout: dict
    :param key: metric value title
    :type key: str
    :param value: metric value
    :type value: anything but str (excluding string representations of bools/Nones)

    :return: Dict containing count of a given metric, value to check if count was made
    :rtype: dict, bool
    """
    if value in ("true", "false", "True", "False", True, False) and isinstance(value, (bool, str)):
        if value in ("true", "True", True):
            qc_layout[f"{key}_counts"] = 1
        else:
            qc_layout[f"{key}_counts"] = 0
        offended = False
    return qc_layout, offended


def count_none(qc_layout, key, value, offended):
    """
    Maps None values to a count of 0.

    :param qc_layout: Dict containing count of a given metric
    :type qc_layout: dict
    :param key: metric value title
    :type key: str
    :param value: metric value
    :type value: anything but str (excluding string representations of bools/Nones)

    :return: Dict containing count of a given metric, value to check if count was made
    :rtype: dict, bool
    """
    if value in ("None", None):
        qc_layout[f"{key}_counts"] = 0
        offended = False
    return qc_layout, offended


def count_datetime(qc_layout, key, dti, offended):
    """
    Maps datetime values to a count of 1.

    :param qc_layout: Dict containing count of a given metric
    :type qc_layout: dict
    :param key: metric value title
    :type key: str
    :param dti: flag representing if datetime was chosen
    :type dti: bool

    :return: Dict containing count of a given metric, value to check if count was made
    :rtype: dict, bool
    """
    if dti:
        qc_layout[f"{key}_counts"] = 1
        offended = False
    return qc_layout, offended


def count_list(qc_layout, key, value, offended):
    """
    Counts list values depending on type (empty, boolean, other)

    :param qc_layout: Dict containing count of a given metric
    :type qc_layout: dict
    :param key: metric value title
    :type key: str
    :param value: metric value
    :type value: anything but str (excluding string representations of bools/Nones)

    :return: Dict containing count of a given metric, value to check if count was made
    :rtype: dict, bool
    """
    if isinstance(value, list):
        # for masks could either be times or bool arrays
        if len(value) == 0:
            counts = 0
        elif not isinstance(value[0], (dict, str)) and set(value) == {
            0,
            1,
        }:  # considers bool case, check for list of ints other than 1 or 0?
            counts = sum(value)
        else:
            counts = len(value)
        qc_layout[f"{key}_counts"] = counts
        offended = False
    return qc_layout, offended


def count_dict(qc_layout, key, value, offended):
    """
    Maps dictionary value to a count of 1.

    :param qc_layout: Dict containing count of a given metric
    :type qc_layout: dict
    :param key: metric value title
    :type key: str
    :param value: metric value
    :type value: anything but str (excluding string representations of bools/Nones)

    :return: Dict containing count of a given metric, value to check if count was made
    :rtype: dict, bool
    """
    if isinstance(value, dict):
        qc_layout[f"{key}_counts"] = 1
        offended = False
    return qc_layout, offended


def count_number(qc_layout, key, value, offended):
    """
    Maps int or float value to a count of the same value.

    :param qc_layout: Dict containing count of a given metric
    :type qc_layout: dict
    :param key: metric value title
    :type key: str
    :param value: metric value
    :type value: anything but str (excluding string representations of bools/Nones)

    :return: Dict containing count of a given metric, value to check if count was made
    :rtype: dict, bool
    """
    # Double check bool since True/False will be instances of (int, float)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        qc_layout[f"{key}_counts"] = value
        offended = False
    return qc_layout, offended


def is_datetime(value):
    """
    Checks to see if value is a datetime string.

    :param value: metric value to be assessed for datetime type
    :type value: any

    :return: Flag to determine if variable is datetime string
    :rtype: bool
    """
    try:
        if not isinstance(value, (int, float)):
            is_dt = UTCDateTime(value)
            dti = isinstance(is_dt, UTCDateTime)
        else:
            dti = False
    except TypeError:
        dti = False
    return dti


def standalone(key, value):
    """
    Records value for variables of interest which are not count data.

    :param key: metric value title
    :type key: str
    :param value: metric value
    :type value: any

    :return: Dict containing count of a given metric
    :rtype: dict
    """
    qc_layout = {f"{key}_value": 0}
    if value is None:
        qc_layout[f"{key}_value"] = 0
    else:
        qc_layout[f"{key}_value"] = value
    return qc_layout


def summary_table_check_exist(db, session_name, network, station, channel, dt_range):
    """
    Checks if the summary table exists. If it does, return None. If it doesn't, return a query to generate it.

    :param db: Database
    :type db: Database
    :param session_name: name of database session
    :type session_name: str
    :param network: network for waveform data
    :type network: str
    :param station: station for waveform data
    :type station: str
    :param channel: channel for waveform data
    :type channel: str
    :param dt_range: datetime range for collected waveform data
    :type dt_range: list of datetime strings

    :return: Query string to generate table if table does not exist
    :rtype: str
    """
    qc_insert = None
    querycheck = """SELECT * FROM summaryReportCountsAndValues 
                    WHERE range = '{dt_range}' AND session = '{sn}' 
                    AND network = '{net}' AND station = '{sta}' 
                    AND channel = '{chan}'""".format(
        dt_range=dt_range, sn=session_name, net=network, sta=station, chan=channel
    )
    db.execute(querycheck)
    selectall_summaryquery = db.fetchall()

    if not selectall_summaryquery:
        qc_insert = """INSERT INTO summaryReportCountsAndValues (created, range, session, 
                        network, station, channel) VALUES ('{created}', '{dt_range}', 
                        '{sn}', '{net}', '{sta}', '{chan}')""".format(
            created=UTCDateTime(time.time()).isoformat(),
            dt_range=dt_range,
            sn=session_name,
            net=network,
            sta=station,
            chan=channel,
        )
    return qc_insert


def uniform_time(starttime, endtime):
    """
    Converts any datetime string to obspy's UTCDateTime type.

    :param starttime: start time of metric data
    :type starttime: datetime str
    :param endtime: end time of metric data
    :type endtime: datetime str

    :return: list of obspy UTCDateTime string or None
    :rtype: list of obspy UTCDatetime string or None
    """
    # Converting to UTCDateTime using the obspy module will use ensure consistent formatting
    # Convert back to string to make human readable
    if starttime and endtime:
        return [UTCDateTime(starttime), UTCDateTime(endtime)]
    else:
        return None


def table_key_qc(key, metric_name):
    """
    Converts given keys to new names. This is used to create unique names
    data when saved to the summary report

    :param key: metric value title
    :type key: str
    :param metric_name: metric name
    :type metric_name: str
    """
    # These values are tied to the summaryReportCountsAndValues table
    # Any changes which happen there, will need to be updated here
    if key == "masks" and metric_name == "deadChanMeanMetric":
        key = "dead_chan_mean_metric_masks"
    if key == "masks" and metric_name == "deadChanADFMetric":
        key = "dead_chan_ADF_metric_masks"
    if key == "masks" and metric_name == "snrMetric":
        key = "snr_masks"
    if key == "count" and metric_name == "repeatedAmplitudeMetric":
        key = "repAmp"
    return key


def pull_root_key(flag):
    if flag.startswith("data_quality") and flag.endswith("counts"):
        root_key = "data_quality"
    elif flag.startswith("io_clock") and flag.endswith("counts"):
        root_key = "io_clock"
    elif flag.startswith("activity") and flag.endswith("counts"):
        root_key = "activity"
    else:
        root_key = None
    return root_key