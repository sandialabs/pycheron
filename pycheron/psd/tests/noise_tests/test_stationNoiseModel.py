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

import pytest

from pycheron.psd.noise.stationNoiseModel import stationNoiseModel


@pytest.fixture(scope="module")
def snm_results_all_types(file_assets):
    stream = file_assets["test_streams"][4]
    results = {}
    types = [
        "envelope10_90",
        "envelope05_95",
        "single10",
        "single90",
        "single05",
        "single95",
    ]
    for type in types:
        results[type] = stationNoiseModel(stream, type=type, plot=False)
    return results


@pytest.mark.parametrize("type", ["envelope10_90", "envelope05_95"])
def test_stationNoiseModel_Has_Correct_Keys_envelope(snm_results_all_types, type):
    result = snm_results_all_types[type]
    result_keys = list(result.keys())
    result_keys.sort()

    expected_keys = [
        "z_station_noiseModel_low",
        "e_station_noiseModel_low",
        "z_station_noiseModel_high",
        "network",
        "metric_name",
        "n_station_noiseModel_low",
        "station",
        "location",
        "e_station_noiseModel_high",
        "n_station_noiseModel_high",
        "type",
    ]
    expected_keys.sort()

    assert result_keys == expected_keys


@pytest.mark.parametrize("type", ["single10", "single90", "single05", "single95"])
def test_stationNoiseModel_Has_Correct_Keys_single(snm_results_all_types, type):
    result = snm_results_all_types[type]
    result_keys = list(result.keys())
    result_keys.sort()

    expected_keys = [
        "station",
        "metric_name",
        "z_station_noiseModel",
        "e_station_noiseModel",
        "n_station_noiseModel",
        "type",
        "location",
        "network",
    ]
    expected_keys.sort()

    assert result_keys == expected_keys


@pytest.mark.parametrize(
    "type",
    [
        "envelope10_90",
        "envelope05_95",
        "single10",
        "single90",
        "single05",
        "single95",
    ],
)
def test_stationNoiseModel_returns_empty_list_for_channels_not_in_trace(snm_results_all_types, type):
    # Trace from Above only has Z channel.
    result = snm_results_all_types[type]
    try:
        assert result["n_station_noiseModel_high"].tolist() == []
        assert result["n_station_noiseModel_low"].tolist() == []
        assert result["e_station_noiseModel_high"].tolist() == []
        assert result["e_station_noiseModel_low"].tolist() == []
    except KeyError:
        # single types of noise models will only have one noiseModel property,
        # not high and low.
        assert result["n_station_noiseModel"].tolist() == []
        assert result["e_station_noiseModel"].tolist() == []


@pytest.mark.parametrize(
    "type",
    [
        "envelope10_90",
        "envelope05_95",
        "single10",
        "single90",
        "single05",
        "single95",
    ],
)
def test_stationNoiseModel_returns_non_empty_list_for_channels_not_in_trace(snm_results_all_types, type):
    result = snm_results_all_types[type]
    assert result["station"] == "ANMO"
    assert result["location"] == "00"
    assert result["network"] == "IU"
    assert result["type"] == type


def test_stationNoiseModel_returns_correct_vals_for_noise_envelope10_90(
    snm_results_all_types,
):
    # Trace only has z channel
    result = snm_results_all_types["envelope10_90"]

    assert result["z_station_noiseModel_high"].tolist()[0] == pytest.approx(-117.490234553, 0.01)
    assert result["z_station_noiseModel_high"].tolist()[-1] == pytest.approx(-74.1643614277, 0.01)

    assert result["z_station_noiseModel_low"].tolist()[0] == pytest.approx(-154.315109347, 0.01)
    assert result["z_station_noiseModel_low"].tolist()[-1] == pytest.approx(-98.5034601199, 0.01)

    assert len(result["z_station_noiseModel_low"]) == 72 == len(result["z_station_noiseModel_high"])


def test_stationNoiseModel_returns_correct_vals_for_noise_envelope05_95(
    snm_results_all_types,
):
    # Trace only has z channel
    result = snm_results_all_types["envelope05_95"]

    assert result["z_station_noiseModel_high"].tolist()[0] == pytest.approx(-117.394009114, 0.01)
    assert result["z_station_noiseModel_high"].tolist()[-1] == pytest.approx(-72.9026636848, 0.01)

    assert result["z_station_noiseModel_low"].tolist()[0] == pytest.approx(-158.821993257, 0.01)
    assert result["z_station_noiseModel_low"].tolist()[-1] == pytest.approx(-100.284149714, 0.01)

    assert len(result["z_station_noiseModel_low"]) == 72 == len(result["z_station_noiseModel_high"])


def test_stationNoiseModel_returns_correct_vals_for_noise_single05(
    snm_results_all_types,
):
    # Trace only has z channel
    result = snm_results_all_types["single05"]
    assert result["z_station_noiseModel"].tolist()[0] == pytest.approx(-158.821993257, 0.01)
    assert result["z_station_noiseModel"].tolist()[-1] == pytest.approx(-100.284149714, 0.01)
    assert len(result["z_station_noiseModel"]) == 72


def test_stationNoiseModel_returns_correct_vals_for_noise_single10(
    snm_results_all_types,
):
    # Trace only has z channel
    result = snm_results_all_types["single10"]
    assert result["z_station_noiseModel"].tolist()[0] == pytest.approx(-154.315109347, 0.01)
    assert result["z_station_noiseModel"].tolist()[-1] == pytest.approx(-98.5034601199, 0.01)
    assert len(result["z_station_noiseModel"]) == 72


def test_stationNoiseModel_returns_correct_vals_for_noise_single90(
    snm_results_all_types,
):
    # Trace only has z channel
    result = snm_results_all_types["single90"]
    assert result["z_station_noiseModel"].tolist()[0] == pytest.approx(-117.490234553, 0.01)
    assert result["z_station_noiseModel"].tolist()[-1] == pytest.approx(-74.1643614277, 0.01)
    assert len(result["z_station_noiseModel"]) == 72


def test_stationNoiseModel_returns_correct_vals_for_noise_single95(
    snm_results_all_types,
):
    # Trace only has z channel
    result = snm_results_all_types["single95"]
    assert result["z_station_noiseModel"].tolist()[0] == pytest.approx(-117.394009114, 0.01)
    assert result["z_station_noiseModel"].tolist()[-1] == pytest.approx(-72.9026636848, 0.01)
    assert len(result["z_station_noiseModel"]) == 72


def test_stationNoiseModel_does_not_error_on_plotting(file_assets):
    stream = file_assets["test_streams"][4]
    stationNoiseModel(stream, type="envelope10_90", plot=True)
