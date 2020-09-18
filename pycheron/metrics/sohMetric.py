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

from future.builtins import *  # NOQA

__all__ = ["sohMetric"]

import obspy.io.mseed.util as util
import obspy
from pycheron.util.logger import Logger
import os
import warnings


def sohMetric(
    st, data_quality=False, activity=False, io_clock=False, logger=None, database=None
):
    # fmt: off
    """
    Extracts accumulated miniSEED quality flags and a measure of timing quality associated with the incoming seismic
    signal. This method will count all set flag bits (data quality, activity, and io/clock) in the fixed section of the
    data header in a Mini-SEED file and returns the total count for each flag type. By default, all flags are returned.
    To return just one set of bit flags (data quality/activity/io_clock), set the optional parameter to True.

    :param st: stream object
    :type st: obspy.core.stream.Stream
    :param data_quality: Return data quality flags
    :type data_quality: bool
    :param activity: Return activity flags
    :type activity: bool
    :param io_clock: Return IO/Clack flags
    :type io_clock: bool
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database

    :return: dictionary with the following keys and types

                * metric_name (`str`)
                * start_time (`str`)
                * end_time (`str`)
                * network (`str`)
                * station (`str`)
                * location (`str`)
                * data_quality_flags (`dict`)
                * activity_flags (`dict`)
                * io_clock_flags (`dict`)

    :rtype: dict

    The miniSEED flags and timing_qual values are described in the SEED manual. [#]_

    Each Stream object contains "accumulators" with counts of the number of times each bit flag was set during the
    parsing of a miniSEED file. Metrics are reported for a subset of these flags

    This method will count all set data quality flag bits in the fixed section of the data header in a Mini-SEED file
    and returns the total count for each flag type.

    **Bit Descriptions**

    .. code-block:: console

        Data Quality
        [Bit 0]	Amplifier saturation detected (station dependent)
        [Bit 1]	Digitizer clipping detected
        [Bit 2]	Spikes detected
        [Bit 3]	Glitches detected
        [Bit 4]	Missing/padded data present
        [Bit 5]	Telemetry synchronization error
        [Bit 6]	A digital filter may be charging
        [Bit 7]	Time tag is questionable

        Activity
        [Bit 0] 	Calibration signals present
        [Bit 1] 	Time correction applied
        [Bit 2] 	Beginning of an event, station trigger
        [Bit 3] 	End of the event, station detriggers
        [Bit 4] 	A positive leap second happened during this record
        [Bit 5] 	A negative leap second happened during this record
        [Bit 6] 	Event in progress

        I/O and Clock
        [Bit 0] 	Station volume parity error possibly present
        [Bit 1] 	Long record read (possibly no problem)
        [Bit 2] 	Short record read (record padded)
        [Bit 3] 	Start of time series
        [Bit 4] 	End of time series
        [Bit 5] 	Clock locked

    **Example**

    .. code-block:: python

        # data with no flags set
        no_flags = 'test/test_data/7a_cabn_bhe.884965.tar.mseed'
        # data with flags set
        flags = 'test/test_data/qualityflags.mseed'

        sohMetric(no_flags) #this will raise an error because there are no data flags associated with it
        >>> "sohMetric: No SOH flags associated with 7A.CABN"
        sohMetric(flags) #can set optional parameters to True to return only certain flag groups
        >>> {'io_clock_flags': {'start_time_series': 0, 'end_time_series': 0, 'long_record_read': 0, 'clock_locked': 0, 'short_record_read': 0, 'station_volume': 0}, 'data_quality_flags': {'digital_filter_charging': 3, 'spikes': 7, 'suspect_time_tag': 2, 'missing_padded_data': 5, 'digitizer_clipping': 8, 'glitches': 6, 'amplifier_saturation': 9, 'telemetry_sync_error': 4}, 'metric_name': 'sohMetric', 'activity_flags': {'event_begin': 0, 'negative_leap': 0, 'event_in_progress': 0, 'event_end': 0, 'time_correction_applied': 0, 'calibration_signal': 0, 'positive_leap': 0}}

    .. rubric:: Footnotes

    .. [#] http://www.fdsn.org/seed_manual/SEEDManual_V2.4.pdf

    """
    # fmt: on
    warnings.filterwarnings("ignore")

    if logger == None:
        logger = Logger(None)

    st_filename = os.path.dirname(__file__) + "/temp.mseed"
    st.write(st_filename)

    try:
        flags = obspy.io.mseed.util.get_flags(st_filename)

        dq = {}

        for k, v in list(flags["data_quality_flags_counts"].items()):
            dq.update({str(k): v})

        ac = {}
        for k, v in list(flags["activity_flags_counts"].items()):
            ac.update({str(k): v})

        io = {}
        for k, v in list(flags["io_and_clock_flags_counts"].items()):
            io.update({str(k): v})

        if data_quality == True:
            return dq
        if activity == True:
            return ac
        if io_clock == True:
            return io

        all = {
            "metric_name": "sohMetric",
            "network": st[0].stats.network,
            "station": st[0].stats.station,
            "location": st[0].stats.location,
            "start_time": st[0].stats.starttime.isoformat(),
            "end_time": st[0].stats.endtime.isoformat(),
        }
        all.update({"data_quality_flags": dq})
        all.update({"activity_flags": ac})
        all.update({"io_clock_flags": io})

        os.remove(st_filename)
        if database is not None:
            database.insert_metric(all)
        return all

    except (TypeError, ValueError, UnboundLocalError):
        st = obspy.read(st_filename)
        sncql = st[0].get_id().split(".")
        logger.error(
            "sohMetric(): No SOH flags associated with " + sncql[0] + "." + sncql[1]
        )
        os.remove(st_filename)
