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

__all__ = ["WfdiscDB", "SQLAlchemyQueryTool", "generate_conn"]


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
        self.site_colnames = [
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
        self.timebox = {"start": 0, "end": 0}
        self.table_names = {
            "wfdisc": None,
            "site": None,
            "sensor": None,
            "instrument": None,
        }

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

    def read_wfdisc_table_time_bound(self, table_name, date_range):
        """Reads rows which fall between date_range in a wfdisc table
        :param table_name: name of the wfdisc table to be accessed
        :type table_name: str
        :param date_range: range of dates to grabe wfdisc rows from
        :type date_range: dict (where keys are 'start' and 'end', and values are obspy.UTCDateTime)

        :return: Returns query result from query tool
        :rtype: `iterable`
        """
        try:
            result = self.query_tool.time_bound_query(table_name, date_range)
        except Exception as e:
            msg = "Failed to read '.w' files: {e}".format(e=e)
            raise ValueError(msg)
        return result

        # Check if returned table is wfdisc table
        if self.wfdisc_colnames != list(result.keys()):
            msg = "Table returned was not a wfdisc table. Expected {wfdisc_colnames}, but got {other_colnames}".format(
                wfdisc_colnames=self.wfdisc_colnames, other_colnames=list(result.keys())
            )
            raise ValueError(msg)
        return result

    def read_site_table_location(self, table_name, sta, timebox):
        """Reads rows which fall between date_range in a wfdisc table
        :param table_name: name of the wfdisc table to be accessed
        :type table_name: str
        :param sta: station to retrieve lat/lon for
        :type sta: dict (where keys are 'start' and 'end', and values are obspy.UTCDateTime)

        :return: Returns query result from query tool
        :rtype: `iterable`
        """
        try:
            result = self.query_tool.location_from_site_table(table_name, sta, timebox)
        except Exception as e:
            msg = "Failed to site table with given time range: {e}".format(e=e)
            raise ValueError(msg)
        return result

        # Check if returned table is site table
        if self.site_colnames != list(result.keys()):
            msg = "Table returned was not a site table. Expected {site_colnames}, but got {other_colnames}".format(
                site_colnames=self.site_colnames, other_colnames=list(result.keys())
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
            print("Converting wfdisc table to pandas dataframe (this may take a while...)")
            for row in wfdisc_table:
                df = df.append(dict(row), ignore_index=True)
        except Exception as e:
            msg = "Failed to convert wfdisc table to pandas dataframe: {e}".format(e=e)
            raise ValueError(msg)
        return df

    def read_instrument_table_resp(self, wfdisc_table_name, sensor_table_name, instr_table_name, sta, cha, timebox):
        """ """
        try:
            result = self.query_tool.resp_file_from_instrument_table(
                wfdisc_table_name,
                sensor_table_name,
                instr_table_name,
                sta,
                cha,
                timebox,
            )
        except Exception as e:
            msg = "Failed to retrieve RESP file location from database: {e}".format(e=e)
            raise ValueError(msg)
        return result

    def set_timebox(self, timebox):
        # TODO T: validate timebox is correct format
        self.timebox = timebox

    def set_table_names(self, table_names):
        # TODO T: Error checking, assure UTCDateTime and keys for variable exist
        self.table_names = table_names


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

    def time_bound_query(self, table_name, date_range):
        # TODO T: remove date_range, fix in testing
        """Crates a timebound query

        :param table_name: string of table name
        :param date_range: dict containing a start and end, where each value is type UTCDateTime
        """
        with self.db_connect.connect() as conn:
            ts = date_range["start"].timestamp
            te = date_range["end"].timestamp

            if date_range["mintime"] == "day":
                mintime = 86400
            elif date_range["mintime"] == "hour":
                mintime = 3600
            elif date_range["mintime"] is None:
                mintime = 0
            else:
                raise ValueError(
                    f"Expected wfdb_mintime from config to be 'day', 'hour', or None. Received: {date_range['mintime']}"
                )

            result = conn.execute(
                f"SELECT * FROM {table_name} WHERE time >= {ts} AND endtime <= {te} and \
                    {table_name}.ENDTIME - {table_name}.TIME >= {mintime}"
            )
        return result

    def location_from_site_table(self, table_name, sta, timebox):
        """
        Retrieves lat/lon info from a site table given some location name
        """
        with self.db_connect.connect() as conn:
            sty, stj = str(timebox["start"].year), str(timebox["start"].julday)
            ety, etj = str(timebox["end"].year), str(timebox["end"].julday)

            # get full julday format expected in site table (year + julday)
            st = sty + stj
            et = ety + etj

            # -1 implies the offdate has not yet happened
            result = conn.execute(
                f"SELECT * FROM {table_name} where sta='{sta}' AND ONDATE <= {st} AND (OFFDATE >= {et} OR OFFDATE = -1)"
            )
        return result

    def resp_file_from_instrument_table(
        self, wfdisc_table_name, sensor_table_name, instr_table_name, sta, cha, timebox
    ):
        """
        Retrieves RESP file from an instrument table given a station and channel
        """
        ts = timebox["start"].timestamp
        te = timebox["end"].timestamp

        with self.db_connect.connect() as conn:
            result = conn.execute(
                f"""
                SELECT
                    {instr_table_name}.DIR,
                    {instr_table_name}.DFILE,
                    {instr_table_name}.RSPTYPE
                FROM {wfdisc_table_name}
                JOIN {sensor_table_name} on
                    {sensor_table_name}.STA='{sta}' and
                    {wfdisc_table_name}.STA ='{sta}' and
                    {sensor_table_name}.CHAN='{cha}' and
                    {wfdisc_table_name}.CHAN='{cha}'
                JOIN {instr_table_name} on
                    {instr_table_name}.INID = {sensor_table_name}.INID
                WHERE {wfdisc_table_name}.TIME >= {ts}
                    and {wfdisc_table_name}.ENDTIME <= {te}
                    and {sensor_table_name}.TIME <= {ts}
                    and {sensor_table_name}.ENDTIME >= {te}
                """
            )
        return result


def generate_conn(orcdb):
    """
    Generates database connection

    :param orcdb: Dictionary of database connection info
    :type orcdb: dict (outlined in example)

    :return: WfdiscDB
    :rtype: WfdiscDB

    Example:
    orcdb = {
             "connect_lib":
                {"cond": cx_Oracle.makedsn, "query_tool": SQLAlchemyQueryTool},
             "connect_dict":
                {"host": "localhost", "port": 1234, "service_name": "myservicename",
                 "user": "myuser", "password": "mypass"},
             "timebox":
                {"start": obspy.UTCDateTime(...), "end": obspy.UTCDateTime(...)},
             "table_names":
                {"wfdisc": "WFDISC_TABLE_NAME", site: "SITE_TABLE_NAME"},
            }
    wfdisc_db = generate_conn(orcdb)
    """
    if orcdb:
        # dereference these to improve some readability
        cld = orcdb["connect_lib"]["cond"]
        qtd = orcdb["connect_lib"]["query_tool"]

        cond = cld(
            orcdb["connect_dict"]["host"],
            orcdb["connect_dict"]["port"],
            service_name=orcdb["connect_dict"]["service_name"],
        )
        qt = qtd(
            "oracle",
            orcdb["connect_dict"]["user"],
            orcdb["connect_dict"]["password"],
            cond,
        )
        wfdisc_db = WfdiscDB(qt)

        wfdisc_db.set_timebox(orcdb["timebox"])
        wfdisc_db.set_table_names(
            {
                "wfdisc": orcdb["table_names"]["wfdisc"],
                "site": orcdb["table_names"]["site"],
                "sensor": orcdb["table_names"]["sensor"],
                "instrument": orcdb["table_names"]["instrument"],
            }
        )
    else:
        wfdisc_db = None

    return wfdisc_db
