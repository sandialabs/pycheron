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

import numpy as np
from obspy.core.utcdatetime import UTCDateTime

from pycheron.psd.psdStatistics import psdStatistics


def test_stats_snlcq_is_correct(BHZ_channel_stream_and_psd_list):
    results = psdStatistics(BHZ_channel_stream_and_psd_list[1])
    for index, tr in enumerate(BHZ_channel_stream_and_psd_list[0]):
        snlcq_tr = "{0}.{1}.{2}.{3}.{4}".format(
            str(tr.stats.network),
            str(tr.stats.station),
            str(tr.stats.location),
            str(tr.stats.channel),
            str(tr.stats.mseed.dataquality),
        )
        assert results[index]["snclq"] == snlcq_tr


def test_stats_metric_name_is_correct(BHZ_channel_stream_and_psd_list):
    results = psdStatistics(BHZ_channel_stream_and_psd_list[1])
    assert results[0]["metric_name"] == "psdStatistics"


def test_start_and_end_times_match(BHZ_channel_stream_and_psd_list):
    results = psdStatistics(BHZ_channel_stream_and_psd_list[1])
    assert UTCDateTime(results[0]["start_time"]) == BHZ_channel_stream_and_psd_list[0][0].stats["starttime"]
    assert UTCDateTime(results[0]["end_time"]) == BHZ_channel_stream_and_psd_list[0][0].stats["endtime"]


def test_stats(BHZ_channel_stream_and_psd_list):
    results = psdStatistics(BHZ_channel_stream_and_psd_list[1])
    noise = results[0]["noise_matrix_noise"]
    assert (results[0]["max"] == np.nanmax(noise, axis=0)).all()
    assert (results[0]["min"] == np.nanmin(noise, axis=0)).all()
    assert (results[0]["mean"] == np.nanmean(noise, axis=0)).all()
    assert (results[0]["median"] == np.nanmedian(noise, axis=0)).all()
    assert (results[0]["percent_95"] == np.nanpercentile(noise, 95, axis=0)).all()
    assert (results[0]["percent_90"] == np.nanpercentile(noise, 90, axis=0)).all()
    assert (results[0]["percent_10"] == np.nanpercentile(noise, 10, axis=0)).all()
    assert (results[0]["percent_5"] == np.nanpercentile(noise, 5, axis=0)).all()


def test_all_other_types(BHZ_channel_stream_and_psd_list):
    results = psdStatistics(BHZ_channel_stream_and_psd_list[1])
    tr_results = results[0]

    assert isinstance(tr_results["noise_matrix_noise"], np.ndarray)
    assert isinstance(tr_results["noise_matrix_noise"][0], np.ndarray)
    assert isinstance(tr_results["noise_matrix_noise"][0][0], np.float64)

    assert isinstance(tr_results["noise_matrix_frequency"], np.ndarray)
    assert isinstance(tr_results["noise_matrix_frequency"][0], np.float64)

    assert isinstance(tr_results["pdf_matrix"], np.ndarray)
    assert isinstance(tr_results["pdf_matrix"][0], np.ndarray)
    assert isinstance(tr_results["pdf_matrix"][0][0], np.float64)

    assert isinstance(tr_results["pdf_bins"], np.ndarray)
    assert isinstance(tr_results["pdf_bins"][0], np.float64)
    # pdf bins are all whole numbers
    assert (np.mod(tr_results["pdf_bins"], 1) == 0).all()

    assert isinstance(tr_results["mode"], np.ndarray)
    assert isinstance(tr_results["mode"][0], np.float64)
    # modes are also all whole numbers
    assert (np.mod(tr_results["mode"], 1) == 0).all()

    assert isinstance(tr_results["nlnm"], np.ndarray)
    assert isinstance(tr_results["nlnm"][0], np.float64)

    assert isinstance(tr_results["nhnm"], np.ndarray)
    assert isinstance(tr_results["nhnm"][0], np.float64)

    assert isinstance(tr_results["percent_below_nlnm"], np.ndarray)
    assert isinstance(tr_results["percent_below_nlnm"][0], np.float64)

    assert isinstance(tr_results["percent_above_nhnm"], np.ndarray)
    assert isinstance(tr_results["percent_above_nhnm"][0], np.float64)


def test_error_handling():
    # These are probably errors and should probably throw a ValueError
    assert psdStatistics(1) == []
    assert psdStatistics(None) == []
    assert psdStatistics([[1]]) == []
    # These may be valid use cases.
    assert psdStatistics([]) == []
    assert psdStatistics([[]]) == []
    # TODO: Error handling in psdStatistics doesn't catch strings passed in
    # ex. psdStatitics('this is an error')
