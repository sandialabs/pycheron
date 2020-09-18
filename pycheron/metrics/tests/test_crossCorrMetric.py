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

from pycheron.metrics.crossCorrMetric import crossCorrMetric
import pycheron.metrics.tests.utils as utils
import pytest


expected_output_lowpass = {
    "peak_correlation": 0.34428011444724643,
    "snclq1_starttime": "2013-11-12T07:09:45",
    "peak_lag": 0.008562450190019734,
    "snclq2_endtime": "2013-11-12T07:19:48",
    "metric_name": "crossCorrMetric",
    "snclq2_starttime": "2013-11-12T07:09:48",
    "snclq": "NM.SLM.00.BHZ.M:NM.SLM.00.BHZ.M",
    "snclq1_endtime": "2013-11-12T07:19:45",
}

expected_output_bandpass = {
    "peak_correlation": -0.06168233613200611,
    "snclq1_starttime": "2013-11-12T07:09:45",
    "peak_lag": -0.0014609637146398536,
    "snclq2_endtime": "2013-11-12T07:19:48",
    "metric_name": "crossCorrMetric",
    "snclq2_starttime": "2013-11-12T07:09:48",
    "snclq": "NM.SLM.00.BHZ.M:NM.SLM.00.BHZ.M",
    "snclq1_endtime": "2013-11-12T07:19:45",
}

expected_output_envelope = {
    "peak_correlation": 0.4353872273981612,
    "snclq1_starttime": "2013-11-12T07:09:45",
    "peak_lag": 0.0108779970188069,
    "snclq2_endtime": "2013-11-12T07:19:48",
    "metric_name": "crossCorrMetric",
    "snclq2_starttime": "2013-11-12T07:09:48",
    "snclq": "NM.SLM.00.BHZ.M:NM.SLM.00.BHZ.M",
    "snclq1_endtime": "2013-11-12T07:19:45",
}

expected_output_int_dec = {
    "peak_correlation": 0.4353872273981612,
    "snclq1_starttime": "2013-11-12T07:09:45",
    "peak_lag": 0.0108779970188069,
    "snclq2_endtime": "2013-11-12T07:19:48",
    "metric_name": "crossCorrMetric",
    "snclq2_starttime": "2013-11-12T07:09:48",
    "snclq": "NM.SLM.00.BHZ.M:NM.SLM.00.BHZ.M",
    "snclq1_endtime": "2013-11-12T07:19:45",
}


expected_output_highpass = {
    "peak_correlation": -0.0069495221096198415,
    "snclq1_starttime": "2013-11-12T07:09:45",
    "peak_lag": 0.0010015403651087618,
    "snclq2_endtime": "2013-11-12T07:19:48",
    "metric_name": "crossCorrMetric",
    "snclq2_starttime": "2013-11-12T07:09:48",
    "snclq": "NM.SLM.00.BHZ.M:NM.SLM.00.BHZ.M",
    "snclq1_endtime": "2013-11-12T07:19:45",
}


expected_output_bandstop = {
    "peak_correlation": 0.4167031487352236,
    "snclq1_starttime": "2013-11-12T07:09:45",
    "peak_lag": 0.010539859005869978,
    "snclq2_endtime": "2013-11-12T07:19:48",
    "metric_name": "crossCorrMetric",
    "snclq2_starttime": "2013-11-12T07:09:48",
    "snclq": "NM.SLM.00.BHZ.M:NM.SLM.00.BHZ.M",
    "snclq1_endtime": "2013-11-12T07:19:45",
}

expected_output_lowpass_cheby_2 = {
    "peak_correlation": 0.9609034963481528,
    "snclq1_starttime": "2013-11-12T07:09:45",
    "peak_lag": 0.024019917586936914,
    "snclq2_endtime": "2013-11-12T07:19:48",
    "metric_name": "crossCorrMetric",
    "snclq2_starttime": "2013-11-12T07:09:48",
    "snclq": "NM.SLM.00.BHZ.M:NM.SLM.00.BHZ.M",
    "snclq1_endtime": "2013-11-12T07:19:45",
}


@pytest.mark.parametrize(
    "filter, expected_output, corners, freq_max, freq_min, max_order",
    [
        ("lowpass", expected_output_lowpass, 2, 100.0, 0.1, None),
        ("bandpass", expected_output_bandpass, 2, 100.0, 0.1, None),
        ("envelope", expected_output_envelope, 2, 100.0, 0.1, None),
        ("integer_decimation", expected_output_int_dec, 2, 100.0, 0.1, None),
        ("highpass", expected_output_highpass, 4, 10.0, 0.1, None),
        ("bandstop", expected_output_bandstop, 4, 10.0, 0.1, None),
        ("bandstop", expected_output_bandstop, 4, 10.0, 0.1, None),
        ("lowpass_cheby_2", expected_output_lowpass_cheby_2, 4, 10.0, 0.1, 2),
    ],
)
def test_crossCorrMetric(
    filter, expected_output, corners, freq_max, freq_min, max_order
):
    st1 = utils.get_stream_from_client(
        "NM", "SLM", "00", "BHZ", "2013-11-12T07:09:45.000", 600
    )
    st2 = utils.get_stream_from_client(
        "NM", "SLM", "00", "BHZ", "2013-11-12T07:09:48", 600
    )
    tr1 = st1[0]
    tr2 = st2[0]
    cross_corr = crossCorrMetric(
        tr1,
        tr2,
        freqmin=freq_min,
        corners=corners,
        filt=filter,
        freqmax=freq_max,
        maxorder=max_order,
    )
    utils.compare_dicts(expected_output, cross_corr)


def test_crossCorrMetricSamplingRate1():
    st1 = utils.get_stream_from_client(
        "NM", "SLM", "00", "BHZ", "2013-11-12T07:09:45.000", 600
    )
    st2 = utils.get_stream_from_client(
        "NM", "SLM", "00", "BHZ", "2013-11-12T07:09:48", 600
    )
    tr1 = st1[0]
    tr2 = st2[0]
    tr1.stats.sampling_rate = 0.9
    cross_corr = crossCorrMetric(tr1, tr2, freqmin=0.1, corners=2)
    assert cross_corr is None


def test_crossCorrMetricSamplingRate2():
    st1 = utils.get_stream_from_client(
        "NM", "SLM", "00", "BHZ", "2013-11-12T07:09:45.000", 600
    )
    st2 = utils.get_stream_from_client(
        "NM", "SLM", "00", "BHZ", "2013-11-12T07:09:48", 600
    )
    tr1 = st1[0]
    tr2 = st2[0]
    tr2.stats.sampling_rate = 0.1
    cross_corr = crossCorrMetric(tr1, tr2, freqmin=0.1, corners=2)
    assert cross_corr is None
