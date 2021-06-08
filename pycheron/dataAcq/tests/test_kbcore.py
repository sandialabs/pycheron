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
from pycheron.dataAcq import kbcore


@pytest.fixture
def site_and_sitechan_assets():
    return {
        "dir": "test_data/",
        "site": "t_test.site",
        "sitechan": "t_test.sitechan",
    }


@pytest.fixture
def site_colnames():
    return [
        "sta",
        "ondate",
        "offdate",
        "lat",
        "lon",
        "elev",
        "staname",
        "statype",
        "refsta",
        "dnorth",
        "deast",
        "lddate",
    ]


@pytest.fixture
def sitechan_colnames():
    return [
        "sta",
        "chan",
        "ondate",
        "chanid",
        "offdate",
        "ctype",
        "edepth",
        "hang",
        "vang",
        "descrip",
        "lddate",
    ]


def test_read_kb_core_site(site_and_sitechan_assets, site_colnames):
    valid_df = kbcore.read_kb_core(site_and_sitechan_assets["dir"] + site_and_sitechan_assets["site"], "site")

    assert valid_df.columns.tolist() == site_colnames
    assert valid_df.empty is False
    assert valid_df.dropna().empty is False


def test_read_kb_core_sitechan(site_and_sitechan_assets, sitechan_colnames):
    valid_df = kbcore.read_kb_core(
        site_and_sitechan_assets["dir"] + site_and_sitechan_assets["sitechan"],
        "sitechan",
    )

    assert valid_df.columns.tolist() == sitechan_colnames
    assert valid_df.empty is False
    assert valid_df.dropna().empty is False


def test_format_kb_core_site(site_and_sitechan_assets, site_colnames):
    with open(site_and_sitechan_assets["dir"] + site_and_sitechan_assets["site"]) as kb:
        for line in kb:
            valid_line = kbcore._format_kb_core_site(line)
            assert valid_line.columns.tolist() == site_colnames
            assert valid_line.empty is False
            assert valid_line.dropna().empty is False


def test_format_kb_core_sitechan(site_and_sitechan_assets, sitechan_colnames):
    with open(site_and_sitechan_assets["dir"] + site_and_sitechan_assets["sitechan"]) as kb:
        for line in kb:
            valid_line = kbcore._format_kb_core_sitechan(line)
            assert valid_line.columns.tolist() == sitechan_colnames
            assert valid_line.empty is False
            assert valid_line.dropna().empty is False


def test_prepare_for_site_or_sitechan(site_colnames, sitechan_colnames):
    valid_site_cols, valid_site_fun = kbcore._prepare_for_site_or_sitechan("site")
    valid_sitechan_cols, valid_sitechan_fun = kbcore._prepare_for_site_or_sitechan("sitechan")

    assert valid_site_cols == site_colnames
    assert callable(valid_site_fun)

    assert valid_sitechan_cols == sitechan_colnames
    assert callable(valid_sitechan_fun)

    with pytest.raises(ValueError):
        invalid_cols, invalid_fun = kbcore._prepare_for_site_or_sitechan("random string")
