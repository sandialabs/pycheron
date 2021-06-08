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
import pytest
import numpy as np
from obspy.clients.fdsn import Client
from pycheron.plotting.dailyPdfPlot import dailyPdfplots
from pycheron.db.sqllite_db import Database
from pycheron.test.create_test_db import create_test_db
from obspy import UTCDateTime
import matplotlib.pyplot as plt

# stand in for enModel and zModel for networkNoiseModel net_dict
test_model = [
    190.2731384004353,
    174.48123722644115,
    159.99999999999991,
    146.72064691274733,
    134.54342644059426,
    123.3768660326352,
    113.13708498984757,
    103.74716437208073,
    95.13656920021765,
    87.24061861322058,
    79.99999999999996,
    73.36032345637366,
    67.27171322029713,
    61.6884330163176,
    56.56854249492378,
    51.873582186040366,
    47.568284600108825,
    43.62030930661029,
    39.99999999999998,
    36.68016172818683,
    33.635856610148565,
    30.8442165081588,
    28.28427124746189,
    25.936791093020183,
    23.784142300054413,
    21.810154653305144,
    19.99999999999999,
    18.340080864093416,
    16.817928305074282,
    15.422108254079408,
    14.142135623730947,
    12.968395546510093,
    11.89207115002721,
    10.905077326652576,
    10.0,
    9.170040432046711,
    8.408964152537145,
    7.711054127039704,
    7.071067811865474,
    6.484197773255047,
    5.946035575013605,
    5.452538663326288,
    5.0,
    4.585020216023356,
    4.204482076268572,
    3.855527063519852,
    3.535533905932737,
    3.2420988866275233,
    2.9730177875068025,
    2.726269331663144,
    2.5,
    2.292510108011678,
    2.102241038134286,
    1.927763531759926,
    1.7677669529663684,
    1.6210494433137617,
    1.4865088937534012,
    1.363134665831572,
    1.25,
    1.146255054005839,
    1.051120519067143,
    0.963881765879963,
    0.8838834764831842,
    0.8105247216568808,
    0.7432544468767006,
    0.681567332915786,
    0.625,
    0.5731275270029195,
    0.5255602595335715,
    0.4819408829399815,
    0.4419417382415921,
    0.4052623608284404,
    0.3716272234383503,
    0.340783666457893,
    0.3125,
    0.28656376350145973,
    0.2627801297667858,
    0.24097044146999075,
    0.22097086912079605,
    0.2026311804142202,
    0.18581361171917515,
    0.1703918332289465,
    0.15625,
    0.14328188175072987,
    0.1313900648833929,
    0.12048522073499537,
    0.11048543456039803,
    0.1013155902071101,
    0.09290680585958758,
    0.08519591661447325,
    0.078125,
    0.07164094087536493,
    0.06569503244169644,
    0.06024261036749766,
    0.05524271728019901,
    0.050657795103555045,
    0.046453402929793775,
    0.04259795830723661,
]


@pytest.fixture
def plotting_assets():
    # set up stream
    # stations/networks/locations/channels may differ, but these values can stay the same
    client = Client("IRIS")
    starttime = UTCDateTime("2012-12-12T00:00:00.000")
    endtime = UTCDateTime("2012-12-16T00:00:00.000")

    # set up database
    create_test_db()
    db_path = "test_data/pycheron_test_db.db"
    db = Database(db_path, "WFpych", True)
    network = "UU"

    # names required for plotting
    f_name_grid = "magrid"
    f_name_line = "maline"

    # User defined networkNoiseModel
    net_dict = {
        "psdStats": {"noise_matrix_frequency": [], "mode": [], "snclq": []},
        "network": "AK",
        "stations_z": ["GHO"],
        "enModel": [test_model, [np.nan] * len(test_model)],
        "metric_name": "networkNoiseModel",
        "stations_en": ["GHO"],
        "zModel": [
            test_model,
            [
                -102.0,
                -103.0,
                -105.0,
                -104.0,
                -104.0,
                -105.0,
                -107.0,
                -107.0,
                -107.0,
                -108.0,
                -108.0,
                -109.0,
                -109.0,
                -110.0,
                -111.0,
                -111.0,
                -112.0,
                -112.0,
                -113.0,
                -113.0,
                -113.0,
                -114.0,
                -114.0,
                -115.0,
                -115.0,
                -116.0,
                -116.0,
                -117.0,
                -117.0,
                -118.0,
                -118.0,
                -119.0,
                -119.0,
                -119.0,
                -119.0,
                -119.0,
                -119.0,
                -117.0,
                -117.0,
                -117.0,
                -117.0,
                -117.0,
                -119.0,
                -119.0,
                -120.0,
                -120.0,
                -121.0,
                -122.0,
                -122.0,
                -122.0,
                -123.0,
                -123.0,
                -123.0,
                -124.0,
                -124.0,
                -124.0,
                -124.0,
                -124.0,
                -125.0,
                -125.0,
                -125.0,
                -125.0,
                -125.0,
                -125.0,
                -125.0,
                -125.0,
                -125.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -126.0,
                -127.0,
                -128.0,
                -128.0,
                -125.0,
                -102.0,
            ],
        ],
    }

    return {
        "client": client,
        "starttime": starttime,
        "endtime": endtime,
        "db": db,
        "db_path": db_path,
        "f_name_grid": f_name_grid,
        "f_name_line": f_name_line,
        "network": network,
        "net_dict": net_dict,
    }


@pytest.mark.parametrize("model", ["gsn", "nlnm", "network"])
def test_db_dailyPdfPlots(plotting_assets, model):

    dailyp = dailyPdfplots(
        plotting_assets["db"],
        model=model,
        f_name_grid=plotting_assets["f_name_grid"],
        f_name_line=plotting_assets["f_name_line"],
        network="UU",
    )
    assert dailyp is not None
    assert os.path.exists(plotting_assets["f_name_grid"] + "_HHE.png")
    assert os.path.exists(plotting_assets["f_name_line"] + "_HHE.png")
    os.remove(plotting_assets["db_path"])
    os.remove(plotting_assets["f_name_grid"] + "_HHE.png")
    os.remove(plotting_assets["f_name_line"] + "_HHE.png")
    plt.close("all")


@pytest.mark.parametrize("model", ["gsn", "nlnm", "network"])
def test_db_database_arg_dailyPdfPlots(plotting_assets, model):
    dailyp_database_arg = dailyPdfplots(
        plotting_assets["db"],
        model=model,
        f_name_grid=plotting_assets["f_name_grid"],
        f_name_line=plotting_assets["f_name_line"],
        network="UU",
        database=plotting_assets["db"],
    )
    assert dailyp_database_arg is not None
    assert os.path.exists(plotting_assets["f_name_grid"] + "_HHE.png")
    assert os.path.exists(plotting_assets["f_name_line"] + "_HHE.png")
    os.remove(plotting_assets["db_path"])
    os.remove(plotting_assets["f_name_grid"] + "_HHE.png")
    os.remove(plotting_assets["f_name_line"] + "_HHE.png")
    plt.close("all")


@pytest.mark.parametrize("model", ["gsn", "nlnm", "network"])
def test_db_net_dict_arg_dailyPdfPlots(plotting_assets, model):
    dailyp_net_dict_arg = dailyPdfplots(
        plotting_assets["db"],
        model=model,
        f_name_grid=plotting_assets["f_name_grid"],
        f_name_line=plotting_assets["f_name_line"],
        network="UU",
        net_dict=plotting_assets["net_dict"],
    )
    assert dailyp_net_dict_arg is not None
    assert os.path.exists(plotting_assets["f_name_grid"] + "_HHE.png")
    assert os.path.exists(plotting_assets["f_name_line"] + "_HHE.png")
    os.remove(plotting_assets["db_path"])
    os.remove(plotting_assets["f_name_grid"] + "_HHE.png")
    os.remove(plotting_assets["f_name_line"] + "_HHE.png")
    plt.close("all")


@pytest.mark.parametrize("model", ["gsn", "nlnm", "network"])
def test_db_database_and_net_dict_args_dailyPdfPlots(plotting_assets, model):
    dailyp_database_and_net_dict_args = dailyPdfplots(
        plotting_assets["db"],
        model=model,
        f_name_grid=plotting_assets["f_name_grid"],
        f_name_line=plotting_assets["f_name_line"],
        network="UU",
        database=plotting_assets["db"],
        net_dict=plotting_assets["net_dict"],
    )
    assert dailyp_database_and_net_dict_args is not None
    assert os.path.exists(plotting_assets["f_name_grid"] + "_HHE.png")
    assert os.path.exists(plotting_assets["f_name_line"] + "_HHE.png")
    os.remove(plotting_assets["db_path"])
    os.remove(plotting_assets["f_name_grid"] + "_HHE.png")
    os.remove(plotting_assets["f_name_line"] + "_HHE.png")
    plt.close("all")


@pytest.mark.parametrize("model", ["gsn", "nlnm", "network"])
def test_st_dailyPdfPlots(plotting_assets, model):
    chan = "BNZ" if model == "network" else "BHE"
    st = plotting_assets["client"].get_waveforms(
        "AK", "GHO", "", chan, plotting_assets["starttime"], plotting_assets["endtime"]
    )

    dailyp = dailyPdfplots(
        st,
        model=model,
        f_name_grid=plotting_assets["f_name_grid"],
        f_name_line=plotting_assets["f_name_line"],
    )
    assert dailyp is not None
    assert os.path.exists(plotting_assets["f_name_grid"] + "_" + chan + ".png")
    assert os.path.exists(plotting_assets["f_name_line"] + "_" + chan + ".png")
    os.remove(plotting_assets["db_path"])
    os.remove(plotting_assets["f_name_grid"] + "_" + chan + ".png")
    os.remove(plotting_assets["f_name_line"] + "_" + chan + ".png")
    plt.close("all")


@pytest.mark.parametrize("model", ["gsn", "nlnm", "network"])
def test_st_database_arg_dailyPdfPlots(plotting_assets, model):
    chan = "BNZ" if model == "network" else "BHE"
    st = plotting_assets["client"].get_waveforms(
        "AK", "GHO", "", chan, plotting_assets["starttime"], plotting_assets["endtime"]
    )

    dailyp_database_arg = dailyPdfplots(
        st,
        model=model,
        f_name_grid=plotting_assets["f_name_grid"],
        f_name_line=plotting_assets["f_name_line"],
        network="AK",
        database=plotting_assets["db"],
    )
    assert dailyp_database_arg is not None
    assert os.path.exists(plotting_assets["f_name_grid"] + "_" + chan + ".png")
    assert os.path.exists(plotting_assets["f_name_line"] + "_" + chan + ".png")
    os.remove(plotting_assets["db_path"])
    os.remove(plotting_assets["f_name_grid"] + "_" + chan + ".png")
    os.remove(plotting_assets["f_name_line"] + "_" + chan + ".png")
    plt.close("all")


@pytest.mark.parametrize("model", ["gsn", "nlnm", "network"])
def test_st_net_dict_arg_dailyPdfPlots(plotting_assets, model):
    chan = "BNZ" if model == "network" else "BHE"
    st = plotting_assets["client"].get_waveforms(
        "AK", "GHO", "", chan, plotting_assets["starttime"], plotting_assets["endtime"]
    )

    dailyp_net_dict_arg = dailyPdfplots(
        st,
        model=model,
        f_name_grid=plotting_assets["f_name_grid"],
        f_name_line=plotting_assets["f_name_line"],
        network="AK",
        net_dict=plotting_assets["net_dict"],
    )
    assert dailyp_net_dict_arg is not None
    assert os.path.exists(plotting_assets["f_name_grid"] + "_" + chan + ".png")
    assert os.path.exists(plotting_assets["f_name_line"] + "_" + chan + ".png")
    os.remove(plotting_assets["db_path"])
    os.remove(plotting_assets["f_name_grid"] + "_" + chan + ".png")
    os.remove(plotting_assets["f_name_line"] + "_" + chan + ".png")
    plt.close("all")


@pytest.mark.parametrize("model", ["gsn", "nlnm", "network"])
def test_st_database_and_net_dict_args_dailyPdfPlots(plotting_assets, model):
    chan = "BNZ" if model == "network" else "BHE"
    st = plotting_assets["client"].get_waveforms(
        "AK", "GHO", "", chan, plotting_assets["starttime"], plotting_assets["endtime"]
    )

    dailyp_database_and_net_dict_args = dailyPdfplots(
        st,
        model=model,
        f_name_grid=plotting_assets["f_name_grid"],
        f_name_line=plotting_assets["f_name_line"],
        network="AK",
        database=plotting_assets["db"],
        net_dict=plotting_assets["net_dict"],
    )
    assert dailyp_database_and_net_dict_args is not None
    assert os.path.exists(plotting_assets["f_name_grid"] + "_" + chan + ".png")
    assert os.path.exists(plotting_assets["f_name_line"] + "_" + chan + ".png")
    os.remove(plotting_assets["db_path"])
    os.remove(plotting_assets["f_name_grid"] + "_" + chan + ".png")
    os.remove(plotting_assets["f_name_line"] + "_" + chan + ".png")
    plt.close("all")
