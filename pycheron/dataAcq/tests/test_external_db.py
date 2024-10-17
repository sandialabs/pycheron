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

import mock
import pytest
import datetime
from pycheron.dataAcq.external_db import WfdiscDB, SQLAlchemyQueryTool
from sqlalchemy.exc import NoSuchModuleError


class MockTable:
    def __init__(self, row):
        self.row = row

    def keys(self):
        return [x[0] for x in self.row[0]]

    def __iter__(self):
        return iter(self.row)


@mock.patch("pycheron.dataAcq.external_db.SQLAlchemyQueryTool")
def test_read_entire_wfdisc_table(mock_sqlalchemy_conn):
    # TODO: mock library becomes unittest.mock in Python 3.5+
    # consider this in the transition from Python 2 -> Python 3
    wfdisc_db = WfdiscDB(mock_sqlalchemy_conn)

    with pytest.raises(ValueError):
        wfdisc_db.read_entire_wfdisc_table("TEST")


@mock.patch("pycheron.dataAcq.external_db.SQLAlchemyQueryTool")
def test_mocked_wfdisc_table_to_dataframe(mock_sqlalchemy_conn):
    mocked_vals = [
        ("sta", "KMBO"),
        ("chan", "BHE"),
        ("time", 1440084800.0),
        ("wfid", 80913458),
        ("chanid", 2023),
        ("jdate", 2015232),
        ("endtime", 1440085859.975),
        ("nsamp", 42400),
        ("samprate", 40.0),
        ("calib", 0.036929),
        ("calper", 1.0),
        ("instype", "STS-2"),
        ("segtype", "o"),
        ("datatype", "e1"),
        ("clip", "-"),
        ("dir", "\\\\mypath\\wf_2015\\232"),
        ("dfile", "KMBO.14_16.2015232.w"),
        ("foff", 133240),
        ("commid", -1),
        ("lddate", datetime.datetime(2016, 1, 14, 16, 5, 9)),
    ]

    mocked_table = MockTable([mocked_vals, mocked_vals, mocked_vals])
    wfdisc_db = WfdiscDB(mock_sqlalchemy_conn)
    df = wfdisc_db.wfdisc_table_to_dataframe(mocked_table)
    assert not df.empty
    assert list(df.columns) == list(mocked_table.keys())


# @mock.patch("pycheron.dataAcq.external_db.SQLAlchemyQueryTool")
# def test_invalid_wfdisc_table_to_dataframe(mock_sqlalchemy_conn):
#     mocked_vals = {""}
#     mocked_table = MockTable([mocked_vals, mocked_vals, mocked_vals])
#     wfdisc_db = WfdiscDB(mock_sqlalchemy_conn)

#     with pytest.raises(ValueError):
#         wfdisc_db.wfdisc_table_to_dataframe(mocked_table)


def test_sql_alchemy_query_tool_creation():
    with pytest.raises(NoSuchModuleError):
        SQLAlchemyQueryTool("test db_type", "test user", "test password", "test host")
