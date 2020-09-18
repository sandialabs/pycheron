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

from pycheron.metrics.sohMetric import sohMetric
import pycheron.metrics.tests.utils as utils


expected_output = {
    "start_time_series": 0,
    "end_time_series": 0,
    "long_record_read": 0,
    "clock_locked": 0,
    "short_record_read": 0,
    "station_volume": 0,
}

expected_output_dq = {
    "digital_filter_charging": 0,
    "spikes": 0,
    "suspect_time_tag": 0,
    "missing_padded_data": 0,
    "digitizer_clipping": 0,
    "glitches": 0,
    "amplifier_saturation": 0,
    "telemetry_sync_error": 0,
}

expected_output_activity = {
    "event_begin": 0,
    "negative_leap": 0,
    "event_in_progress": 0,
    "event_end": 0,
    "time_correction_applied": 0,
    "calibration_signal": 0,
    "positive_leap": 0,
}


expected_output_io_clock = {
    "start_time_series": 0,
    "end_time_series": 0,
    "long_record_read": 0,
    "clock_locked": 0,
    "short_record_read": 0,
    "station_volume": 0,
}


def test_sohMetric():
    st = utils.get_stream_from_file("test_data/qualityflags.mseed")
    soh = sohMetric(st)
    assert len(soh) == 9
    utils.compare_dicts(expected_output, soh["io_clock_flags"])


def test_sohMetricDQ():
    st = utils.get_stream_from_file("test_data/qualityflags.mseed")
    soh = sohMetric(st, data_quality=True)
    assert len(soh) == 8
    utils.compare_dicts(expected_output_dq, soh)


def test_sohMetricActivity():
    st = utils.get_stream_from_file("test_data/qualityflags.mseed")
    soh = sohMetric(st, activity=True)
    assert len(soh) == 7
    utils.compare_dicts(expected_output_activity, soh)


def test_sohMetricIO():
    st = utils.get_stream_from_file("test_data/qualityflags.mseed")
    soh = sohMetric(st, io_clock=True)
    assert len(soh) == 6
    utils.compare_dicts(expected_output_io_clock, soh)
