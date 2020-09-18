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

import pandas as pd
import sqlalchemy


class WfdiscDB:
    """Class responsible for performing queries on a WFDisc table

    **Example**

    .. code-block:: python
        import cx_Oracle

        connect_dict = {
        'host': 'examplehost',
        'port': '5432',
        'sid': 'mysid',
        'user': 'myuser',
        'password': os.environ.get('ORACLE_TEMP_PYCHERON'),
        }

        sid = cx_Oracle.makedsn(connect_dict["host"], connect_dict["port"], sid=connect_dict["sid"])
        orcl_alchemy = SQLAlchemyQueryTool("oracle", connect_dict["user"], connect_dict["password"], sid)
        wfdisc_db = WfdiscDB(orcl_alchemy)
        wfdisc_table = wfdisc_db.read_entire_wfdisc_table('WFDISC_SNL')
        wfdisc_df = wfdisc_db.wfdisc_table_to_dataframe(wfdisc_table)

    """

    def __init__(self, query_tool):
        """
        :param query_tool: Selected library to perform queries.
        :type query_tool: Class (with db_connect variable and table_query method)
        :param wfdisc_colnames: Expected column names for a wfdisc table
        :type wfdisc_colnames: list
        """
        self.query_tool = query_tool
        self.wfdisc_colnames = [
            "sta",
            "chan",
            "time",
            "wfid",
            "chanid",
            "jdate",
            "endtime",
            "nsamp",
            "samprate",
            "calib",
            "calper",
            "instype",
            "segtype",
            "datatype",
            "clip",
            "dir",
            "dfile",
            "foff",
            "commid",
            "lddate",
        ]

    def read_entire_wfdisc_table(self, table_name):
        """Reads all rows in a wfdisc table
        :param table_name: name of the wfdisc table to be accessed
        :type table_name: str

        :return: Returns query result from query tool
        :rtype: `iterable`
        """
        result = self.query_tool.table_query(table_name)

        # Check if returned table is wfdisc table
        if self.wfdisc_colnames != list(result.keys()):
            msg = "Table returned was not a wfdisc table. Expected {wfdisc_colnames}, but got {other_colnames}".format(
                wfdisc_colnames=self.wfdisc_colnames, other_colnames=list(result.keys())
            )
            raise ValueError(msg)
        return result

    def wfdisc_table_to_dataframe(self, wfdisc_table):
        """Converts wfdisc table into a pandas dataframe
        :param wfdisc_table: wfdisc table returned from query tool
        :type wfdisc_table: `iterable`

        :return: Returns a pandas dataframe representation of the wfdisc iterable
        :type: `pandas.DataFrame`
        """
        df = pd.DataFrame(columns=list(wfdisc_table.keys()))
        try:
            print(
                "Converting wfdisc table to pandas dataframe (this may take a while...)"
            )
            for row in wfdisc_table:
                df = df.append(dict(row), ignore_index=True)
        except Exception as e:
            msg = "Failed to convert wfdisc table to pandas dataframe: {e}".format(e=e)
            raise ValueError(msg)
        return df


class SQLAlchemyQueryTool:
    """SQLAlchemy specific Querytool to be used for WfdiscDB.
    All future query tools MUST HAVE THE SAME CLASS VARIABLES AND METHODS!
    """

    def __init__(self, db_type, user, password, host):
        """
        :param db_type: type of database to be accessed (oracle, postgres, etc.)
        :type db_type: str
        :param user: username for database access
        :type user: str
        :param password: password for database access
        :type password: str
        :param host: entire hostname database exists on
        :type host: str
        """
        cstr = sqlalchemy.engine.url.URL(db_type, user, password, host)
        self.db_connect = sqlalchemy.create_engine(
            cstr, convert_unicode=False, pool_recycle=10, pool_size=50, echo=True
        )

    def table_query(self, table_name):
        """Reads all rows in a wfdisc table
        :param table_name: name of the wfdisc table to be accessed
        :type table_name: str

        :return: Returns query result from query tool
        :rtype: `iterable`
        """
        with self.db_connect.connect() as conn:
            table_obj = sqlalchemy.Table(
                table_name,
                sqlalchemy.MetaData(),
                autoload=True,
                autoload_with=self.db_connect,
            )
            s = sqlalchemy.sql.select([table_obj])
            result = conn.execute(s)
        return result
