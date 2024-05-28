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

__all__ = ["dbIntegrityCheck"]

from sqlalchemy import create_engine
from sqlalchemy.engine.url import URL
from sqlalchemy.sql import text
from cx_Oracle import makedsn
from pycheron.util.logger import Logger
from pycheron.db.sqllite_db import Database


def dbIntegrityCheck(
    dialect="oracle",
    user=None,
    pswd=None,
    host=None,
    port=None,
    service_name=None,
    CSSType="CSS",
    table_names=None,
    logger=None,
    database_config=None,
):
    """

    Basic Oracle CSS3.0/KB COre database integrity checks on affiliation, instrument, sensor, site, sitechan, and
    wfdisc tables

    :param dialect: name of the SQLAlchemy dialect, such as sqlite, mysql, postgresql, oracle, or mssql.
                             Code only currently works with oracle (default)
    :type dialect: str
    :param user: user name of the database (default = None)
    :type user: str
    :param pswd: password for the database (default = None)
    :type pswd: str
    :param host: host to connect to (default = None)
    :type host: str
    :param port: port to connect to (default = None)
    :type port: str
    :param service_name: Oracle db service name (default = None)
    :type service_name: str
    :param CSSType: string delineating schema type. Current options or 'CSS', or 'KB Core' (default = 'CSS')
    :type CSSType: str
    :param table_names: names of wfdisc tables being tested for integrity (refer to callPycheronConfig table_names)
    :type table_names: dict
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger

    Input format should provide info to create this information:
        dialect[+driver]://user:password@host:port/service name.

    ** Note **
    Depending on where cx_Oracle client installed, users may have trouble with cx_Oracle driver if it can't find the lib
    dir. May need to initialize cx_Oracle with given path, i.e., cx_Oracle.init_oracle_client(lib_dir="path_to_lib_dir")

    return: basic_db_check: dictionary of metric name and integrity issue counts
    rtype: dict

    Dictionary output is the following:
    'metric_name': Name of metric
    'integrity_results': {
    'Affiliation_EndTime_Integrity_Issues': Number of entries in affiliation table where endtime > 9999999999.999
    'Affiilation_Time_Integrity_Issues': Number of entries in affiliation table where time < -9999999999.999
    'Instrument_Digital_Integrity_Issues': Number of entries in instrument table where digital NOT IN ('d','a','-')
    'Instrument_Inid_Integrity_Issues': Number of entries in instrument table where inid <= 0
    'Instrument_NCalib_Integrity_Issues': Number of entries in instrument table where ncalib = 0
    'Instrument_Ncalper_Integrity_Issues': Number of entries in instrument table where ncalper < 0
    'Instrument_Rsptype_Integrity_Issues': Number of entries in instrument table where rsptype NOT IN
                                           ('paz','fir','fap', 'pazfir', 'evresp', 'sacpzf', '-')
    'Instrument_SampRate_Integrity_Issues': Number of entries in instrument table where samprate <= 0
    'Instrument_Band_Integrity_Issues': Number of entries in instrument table where band NOT IN
                                        ('s','m','i','l','b','h','v','e','r','w','a','f','c','g','-')")
    'Sensor_Calper_Integrity_Issues': Number of entries in sensor table where calper <= 0
    'Sensor_CalRatio_Integrity_Issues': Number of entries in sensor table where calratio <= 0
    'Sensor_Chanid_Integrity_Issues': Number of entries in sensor table where chanid < 0 and chanid != -1
    'Sensor_EndTime_Integrity_Issues': Number of entries in sensor table where endtime > 9999999999.999
    'Sensor_Inid_Integrity_Issues': Number of entries in sensor table where inid < 0 and inid != -1
    'Sensor_Instant_Integrity_Issues': Number of entries in sensor table where instant NOT IN ('y', 'n')
    'Sensor_Jdate_Integrity_Issues': Number of entries in sensor table where jdate < 1901348 or jdate > 3001000
    'Sensor_Time_Integrity_Issues': Number of entries in sensor table where time < -9999999999.999
    'Sensor_Tshift_Integrity_Issues': Number of entries in sensor table where tshift < -9999999999.999
    'Site_DEast_Integrity_Issues': Number of entries in site table where deast < -20000 or deast > 20000
    'Site_DNorth_Integrity_Issues': Number of entries in site table where dnorth < -20000 or dnorth > 20000
    'Site_Elev_Integrity_Issues': Number of entries in site table where elev < -10 or elev > 10
    'Site_Lat_Integrity_Issues': Number of entries in site table where lat < -90 or lat > 90
    'Site_Lon_Integrity_Issues': Number of entries in site table where lon < -180 or lon > 180
    'Site_OffDate_Integrity_Issues': Number of entries in site table where offdate < 1901348 or offdate > 3001000
    'Site_OnDate_Integrity_Issues': Number of entries in site table where ondate < 1800001 or ondate > 3001000
    'Site_StaType_Integrity_Issues': Number of entries in site table where statype NOT IN ('ss', 'ar', '-')
    'Sitechan_Chanid_Integrity_Issues': Number of entries in sitechan table where chanid < 0 and chanid != -1
    'Sitechan_Ctype_Integrity_Issues': Number of entries in sitechan table where ctype NOT IN ('n', 'b', 'i', '-')
    'Sitechan_Hang_Integrity_Issues': Number of entries in sitechan table where hang < 0 or hang > 360
    'Sitechan_Vang_Integrity_Issues': Number of entries in sitechan table where vang < 0 or vang > 90
    'Sitechan_OffDate_Integrity_Issues': Number of entries in sitechan table where offdate < 1901348 or offdate >3001000
    'Sitechan_OnDate_Integrity_Issues': Number of entries in sitechan table where ondate < 1800001 or ondate > 3001000
    'Sitechan_Edepth_Integrity_Issues': Number of entries in sitechan table where edepth < -10
    'Wfdisc_Calib_Integrity_Issues': Number of entries in wfdisc table where calib <= 0
    'Wfdisc_Calper_Integrity_Issues': Number of entries in wfdisc table where calper <= 0
    'Wfdisc_Chanid_Integrity_Issues': Number of entries in wfdisc table where chanid < 0 and chanid != -1
    'Wfdisc_Clip_Integrity_Issues': Number of entries in wfdisc table where clip NOT IN ('c', 'n', '-')
    'Wfdisc_Commid_Integrity_Issues': Number of entries in wfdisc table where commid < 0 and commid != -1
    'Wfdisc_EndTime_Integrity_Issues': Number of entries in wfdisc table where endtime > 9999999999.999
    'Wfdisc_Foff_Integrity_Issues': Number of entries in wfdisc table where foff < 0
    'Wfdisc_Jdate_Integrity_Issues': Number of entries in wfdisc table where jdate < 1901348 and jdate > 3001000
    'Wfdisc_Nsamp_Integrity_Issues': Number of entries in wfdisc table where nsamp <= 0
    'Wfdisc_Samprate_Integrity_Issues': Number of entries in wfdisc table where samprate <= 0
    'Wfdisc_SegType_Integrity_Issues': Number of entries in wfdisc table where segtype NOT IN ('o', 'v', 's', 'd', '-')
    'Wfdisc_Time_Integrity_Issues': Number of entries in wfdisc table where time < -9999999999.999
    'Wfdisc_Wfid_Integrity_Issues': Number of entries in wfdisc table where wfid <= 0
    'Wfdisc_Datatype_Integrity_Issues': Number of entries in wfdisc table where datatype NOT IN
                                        ('t4', 'e1', 's4', 's3', 's2', 'g2', 'i4', '-')
    }
    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Check CSSType input as only allow CSS or KB Core currently
    if CSSType not in ("CSS", "KB Core"):
        logger.error("CSSType must be either CSS or KB Core")
        return

    # Set up dbinfo for create engine to connect to the database
    sn = makedsn(host, port, service_name=service_name)
    cstr = URL(dialect, user, pswd, sn)
    
    # Create engine to connect to the database
    db_connect = create_engine(cstr, convert_unicode=False, pool_recycle=10, pool_size=50, echo=True)
    
    # Connect to the database
    connection = db_connect.connect()

    # Initialize variables
    AffEndTimeCount = None
    AffTimeCount = None
    InstDigitalCount = None
    InstInidCount = None
    InstNcalibCount = None
    InstNcalperCount = None
    InstRsptypeCount = None
    InstSampRateCount = None
    InstBandCount = None
    SensorCalperCount = None
    SensorCalRatioCount = None
    SensorChanidCount = None
    SensorEndTimeCount = None
    SensorInidCount = None
    SensorInstantCount = None
    SensorJdateCount = None
    SensorTimeCount = None
    SensorTshiftCount = None
    SiteDEastCount = None
    SiteDNorthCount = None
    SiteElevCount = None
    SiteLatCount = None
    SiteLonCount = None
    SiteOffDateCount = None
    SiteOnDateCount = None
    SiteStaTypeCount = None
    SitechanChanidCount = None
    SitechanCtypeCount = None
    SitechanHangCount = None
    SitechanVangCount = None
    SitechanOffDateCount = None
    SitechanOndateCount = None
    SitechanEdepthCount = None
    WfiscCalibCount = None
    WfdiscCalperCount = None
    WfdiscChanidCount = None
    WfdiscClipCount = None
    WfdiscCommidCount = None
    WfdiscEndTimeCount = None
    WfdiscFoffCount = None
    WfdiscJdateCount = None
    WfdiscNsampCount = None
    WfdiscSamprateCount = None
    WfdiscSegTypeCount = None
    WfdiscTimeCount = None
    WfdiscWfidCount = None
    WfdiscDatatypeCount = None

    # Now that we have a connection, do general database integrity checks on existing tables
    ######### AFFILIATION TABLE ############################################################
    # Set up queries for basic affiliation table checks ane execute query
    if CSSType == "KB Core":
        # Integrity check on endtime
        AffEndTimeQuery = text(
            f"SELECT count(*) FROM {table_names['affiliation']} WHERE endtime > 9999999999.999"
        )
        AffEndTimeResult = connection.execute(AffEndTimeQuery)
        for affET in AffEndTimeResult:
            AffEndTimeCount = affET[0]
        # Integrity check on time
        AffTimeQuery = text(
            f"SELECT count(*) FROM {table_names['affiliation']} WHERE time < -9999999999.999"
        )
        AffTimeResult = connection.execute(AffTimeQuery)
        for affTime in AffTimeResult:
            AffTimeCount = affTime[0]

    ######### INSTRUMENT TABLE #############################################################
    # Set up queries for basic instrument table checks and execute query
    # Integrity check on digital
    InstDigitalQuery = text(
        f"SELECT count(*) FROM {table_names['instrument']} WHERE digital NOT IN ('d','a','-')"
    )
    InstDigitalResult = connection.execute(InstDigitalQuery)
    for InstDig in InstDigitalResult:
        InstDigitalCount = InstDig[0]
    # Integrity check on inid
    InstInidQuery = text(f"SELECT count(*) FROM {table_names['instrument']} WHERE inid <= 0")
    InstInidResult = connection.execute(InstInidQuery)
    for InstInid in InstInidResult:
        InstInidCount = InstInid[0]
    # Integrity check on ncalib
    InstNcalibQuery = text(f"SELECT count(*) FROM {table_names['instrument']} WHERE ncalib = 0")
    InstNcalibResult = connection.execute(InstNcalibQuery)
    for InstNcalib in InstNcalibResult:
        InstNcalibCount = InstNcalib[0]
    # Integrity check on ncalper
    InstNcalperQuery = text(f"SELECT count(*) FROM {table_names['instrument']} WHERE ncalper <= 0")
    InstNcalperResult = connection.execute(InstNcalperQuery)
    for InstNcalper in InstNcalperResult:
        InstNcalperCount = InstNcalper[0]
    # Integrity check on rsptype
    InstRsptypeQuery = text(
        f"SELECT count(*) FROM {table_names['instrument']} WHERE rsptype NOT IN ('paz','fir','fap', 'pazfir', 'evresp', 'sacpzf', '-')"
    )
    InstRsptypeResult = connection.execute(InstRsptypeQuery)
    for InstRsptype in InstRsptypeResult:
        InstRsptypeCount = InstRsptype[0]
    # Integrity check on sample rate
    InstSampRateQuery = text(f"SELECT count(*) FROM {table_names['instrument']} WHERE samprate <= 0")
    InstSampRateResult = connection.execute(InstSampRateQuery)
    for InstSampRate in InstSampRateResult:
        InstSampRateCount = InstSampRate[0]
    # Integrity check on band
    if CSSType == "CSS":
        InstBandQuery = text(
            f"SELECT count(*) FROM {table_names['instrument']} WHERE band NOT IN ('s','m','i','l','b','h','v', '-')"
        )
        InstBandResult = connection.execute(InstBandQuery)
        for InstBand in InstBandResult:
            InstBandCount = InstBand[0]
    else:
        InstBandQuery = text(
            f"SELECT count(*) FROM {table_names['instrument']} WHERE band NOT IN ('s','m','i','l','b','h','v', 'e', 'r', 'u', 'w', '-')"
        )
        InstBandResult = connection.execute(InstBandQuery)
        for InstBand in InstBandResult:
            InstBandCount = InstBand[0]

    ######### SENSOR TABLE #################################################################
    # Set up queries for basic sensor table checks and execute query
    # Integrity check on calper
    SensorCalperQuery = text(f"SELECT count(*) FROM {table_names['sensor']} WHERE calper <= 0")
    SensorCalperResult = connection.execute(SensorCalperQuery)
    for SensorCalper in SensorCalperResult:
        SensorCalperCount = SensorCalper[0]
    # Integrity check on calratio
    SensorCalRatioQuery = text(f"SELECT count(*) FROM {table_names['sensor']} WHERE calratio = 0")
    SensorCalRatioResult = connection.execute(SensorCalRatioQuery)
    for SensorCalRatio in SensorCalRatioResult:
        SensorCalRatioCount = SensorCalRatio[0]
    # Integrity check on chanid
    SensorChanidQuery = text(
        f"SELECT count(*) FROM {table_names['sensor']} WHERE chanid <= 0 and chanid != -1"
    )
    SensorChanidResult = connection.execute(SensorChanidQuery)
    for SensorChanid in SensorChanidResult:
        SensorChanidCount = SensorChanid[0]
    # Integrity check on endtime
    SensorEndTimeQuery = text(
        f"SELECT count(*) FROM {table_names['sensor']} WHERE endtime > 9999999999.999"
    )
    SensorEndTimeResult = connection.execute(SensorEndTimeQuery)
    for SensorEndTime in SensorEndTimeResult:
        SensorEndTimeCount = SensorEndTime[0]
    # Integrity check on inid
    SensorInidQuery = text(f"SELECT count(*) FROM {table_names['sensor']} WHERE inid <= 0 and inid != -1")
    SensorInidResult = connection.execute(SensorInidQuery)
    for SensorInid in SensorInidResult:
        SensorInidCount = SensorInid[0]
    # Integrity check on instant
    SensorInstantQuery = text(
        f"SELECT count(*) FROM {table_names['sensor']} WHERE instant NOT IN ('y', 'n')"
    )
    SensorInstantResult = connection.execute(SensorInstantQuery)
    for SensorInstant in SensorInstantResult:
        SensorInstantCount = SensorInstant[0]
    # Integrity check on jdate
    SensorJdateQuery = text(
        f"SELECT count(*) FROM {table_names['sensor']} WHERE jdate < 1901348 or jdate > 3001000"
    )
    SensorJdateResult = connection.execute(SensorJdateQuery)
    for SensorJdate in SensorJdateResult:
        SensorJdateCount = SensorJdate[0]
    # Integrity check on time
    SensorTimeQuery = text(f"SELECT count(*) FROM {table_names['sensor']} WHERE time < -9999999999.999")
    SensorTimeResult = connection.execute(SensorTimeQuery)
    for SensorTime in SensorTimeResult:
        SensorTimeCount = SensorTime[0]
    # Integrity check on tshift
    # Technically there isn't a range specified for CSS but if below this number there's an obvious issue
    SensorTshiftQuery = text(
        f"SELECT count(*) FROM {table_names['sensor']} WHERE tshift <= -9999999999.999"
    )
    SensorTshiftResult = connection.execute(SensorTshiftQuery)
    for SensorTshift in SensorTshiftResult:
        SensorTshiftCount = SensorTshift[0]

    ######### SITE TABLE ####################################################################
    # Set up queries for basic site table checks and execute query
    # Integrity check on deast
    SiteDEastQuery = text(
        f"SELECT count(*) FROM {table_names['site']} WHERE deast < -20000 or deast > 20000"
    )
    SiteDEastResult = connection.execute(SiteDEastQuery)
    for SiteDEast in SiteDEastResult:
        SiteDEastCount = SiteDEast[0]
    # Integrity check on dnorth
    SiteDNorthQuery = text(
        f"SELECT count(*) FROM {table_names['site']} WHERE dnorth < -20000 or dnorth > 20000"
    )
    SiteDNorthResult = connection.execute(SiteDNorthQuery)
    for SiteDNorth in SiteDNorthResult:
        SiteDNorthCount = SiteDNorth[0]
    # Integrity check on elev
    SiteElevQuery = text(f"SELECT count(*) FROM {table_names['site']} WHERE elev < -10 or elev > 10")
    SiteElevResult = connection.execute(SiteElevQuery)
    for SiteElev in SiteElevResult:
        SiteElevCount = SiteElev[0]
    # Integrity check on lat
    SiteLatQuery = text(f"SELECT count(*) FROM {table_names['site']} WHERE lat < -90 or lat > 90")
    SiteLatResult = connection.execute(SiteLatQuery)
    for SiteLat in SiteLatResult:
        SiteLatCount = SiteLat[0]
    # Integrity check on lon
    SiteLonQuery = text(f"SELECT count(*) FROM {table_names['site']} WHERE lon < -180 or lon > 180")
    SiteLonResult = connection.execute(SiteLonQuery)
    for SiteLon in SiteLonResult:
        SiteLonCount = SiteLon[0]
    # Integrity check on offdate
    SiteOffDateQuery = text(
        f"SELECT count(*) FROM {table_names['site']} WHERE offdate < 1901348 or offdate > 3001000"
    )
    SiteOffDateResult = connection.execute(SiteOffDateQuery)
    for SiteOffDate in SiteOffDateResult:
        SiteOffDateCount = SiteOffDate[0]
    # Integrity check on ondate
    SiteOnDateQuery = text(
        f"SELECT count(*) FROM {table_names['site']} WHERE ondate < 1800001 or ondate > 3001000"
    )
    SiteOnDateResult = connection.execute(SiteOnDateQuery)
    for SiteOnDate in SiteOnDateResult:
        SiteOnDateCount = SiteOnDate[0]
    # Integrity check on statype
    SiteStaTypeQuery = text(
        f"SELECT count(*) FROM {table_names['site']} WHERE statype NOT IN ('ss', 'ar', '-')"
    )
    SiteStaTypeResult = connection.execute(SiteStaTypeQuery)
    for SiteStaType in SiteStaTypeResult:
        SiteStaTypeCount = SiteStaType[0]

    ######### SITECHAN TABLE ################################################################
    # Set up queries for basic sitechan table checks and execute query
    # Integrity check on chanid
    SitechanChanidQuery = text(
        f"SELECT count(*) FROM {table_names['sitechan']} WHERE chanid <= 0 and chanid != -1"
    )
    SitechanChanidResult = connection.execute(SitechanChanidQuery)
    for SitechanChanid in SitechanChanidResult:
        SitechanChanidCount = SitechanChanid[0]
    # Integrity check on ctype
    SitechanCtypeQuery = text(
        f"SELECT count(*) FROM {table_names['sitechan']} WHERE ctype NOT IN ('n', 'b', 'i', '-')"
    )
    SitechanCtypeResult = connection.execute(SitechanCtypeQuery)
    for SitechanCtype in SitechanCtypeResult:
        SitechanCtypeCount = SitechanCtype[0]
    # Integrity check on hang
    SitechanHangQuery = text(
        f"SELECT count(*) FROM {table_names['sitechan']} WHERE hang < 0 or hang > 360"
    )
    SitechanHangResult = connection.execute(SitechanHangQuery)
    for SitechanHang in SitechanHangResult:
        SitechanHangCount = SitechanHang[0]
    # Integrity check on vang
    if CSSType == "CSS":
        SitechanVangQuery = text(
            f"SELECT count(*) FROM {table_names['sitechan']} WHERE vang < 0 or vang > 90"
        )
        SitechanVangResult = connection.execute(SitechanVangQuery)
        for SitechanVang in SitechanVangResult:
            SitechanVangCount = SitechanVang[0]
    else:
        SitechanVangQuery = text(
            f"SELECT count(*) FROM {table_names['sitechan']} WHERE vang < 0 or vang > 180"
        )
        SitechanVangResult = connection.execute(SitechanVangQuery)
        for SitechanVang in SitechanVangResult:
            SitechanVangCount = SitechanVang[0]
    # Integrity check on offdate
    SitechanOffDateQuery = text(
        f"SELECT count(*) FROM {table_names['sitechan']} WHERE offdate < 1901348 or offdate > 3001000"
    )
    SitechanOffDateResult = connection.execute(SitechanOffDateQuery)
    for SitechanOffDate in SitechanOffDateResult:
        SitechanOffDateCount = SitechanOffDate[0]
    # Integrity check on ondate
    SitechanOndateQuery = text(
        f"SELECT count(*) FROM {table_names['sitechan']} WHERE ondate < 1800001 or ondate > 3001000"
    )
    SitechanOndateResult = connection.execute(SitechanOndateQuery)
    for SitechanOndate in SitechanOndateResult:
        SitechanOndateCount = SitechanOndate[0]
    # Integrity check on edepth
    SitechanEdepthQuery = text(f"SELECT count(*) FROM {table_names['sitechan']} WHERE edepth < 0")
    SitechanEdepthResult = connection.execute(SitechanEdepthQuery)
    for SitechanEdepth in SitechanEdepthResult:
        SitechanEdepthCount = SitechanEdepth[0]

    ######### WFDISC TABLE ####################################################################
    # Set up queries for basic wfdisc table checks and execute query
    # Integrity check on calib
    WfdiscCalibQuery = text(f"SELECT count(*) FROM {table_names['wfdisc']} WHERE calib <= 0")
    WfiscCalibResult = connection.execute(WfdiscCalibQuery)
    for WfiscCalib in WfiscCalibResult:
        WfiscCalibCount = WfiscCalib[0]
    # Integrity check on calper
    WfdiscCalperQuery = text(f"SELECT count(*) FROM {table_names['wfdisc']} WHERE calper <= 0")
    WfdiscCalperResult = connection.execute(WfdiscCalperQuery)
    for WfdiscCalper in WfdiscCalperResult:
        WfdiscCalperCount = WfdiscCalper[0]
    # Integrity check on chanid
    WfdiscChanidQuery = text(
        f"SELECT count(*) FROM {table_names['wfdisc']} WHERE chanid <= 0 and chanid != -1"
    )
    WfdiscChanidResult = connection.execute(WfdiscChanidQuery)
    for WfdiscChanid in WfdiscChanidResult:
        WfdiscChanidCount = WfdiscChanid[0]
    # Integrity check on clip
    WfdiscClipQuery = text(
        f"SELECT count(*) FROM {table_names['wfdisc']} WHERE clip NOT IN ('c', 'n', '-')"
    )
    WfdiscClipResult = connection.execute(WfdiscClipQuery)
    for WfdiscClip in WfdiscClipResult:
        WfdiscClipCount = WfdiscClip[0]
    # Integrity check on commid
    WfdiscCommidQuery = text(
        f"SELECT count(*) FROM {table_names['wfdisc']} WHERE commid <= 0 and commid != -1"
    )
    WfdiscCommidResult = connection.execute(WfdiscCommidQuery)
    for WfdiscCommid in WfdiscCommidResult:
        WfdiscCommidCount = WfdiscCommid[0]
    # Integrity check on endtime
    WfdiscEndTimeQuery = text(
        f"SELECT count(*) FROM {table_names['wfdisc']} WHERE endtime > 9999999999.999"
    )
    WfdiscEndTimeResult = connection.execute(WfdiscEndTimeQuery)
    for WfdiscEndTime in WfdiscEndTimeResult:
        WfdiscEndTimeCount = WfdiscEndTime[0]
    # Integrity check on foff
    WfdiscFoffQuery = text(f"SELECT count(*) FROM {table_names['wfdisc']} WHERE foff < 0")
    WfdiscFoffResult = connection.execute(WfdiscFoffQuery)
    for WfdiscFoff in WfdiscFoffResult:
        WfdiscFoffCount = WfdiscFoff[0]
    # Integrity check on jdate
    WfdiscJdateQuery = text(
        f"SELECT count(*) FROM {table_names['wfdisc']} WHERE jdate < 1901348 and jdate > 3001000"
    )
    WfdiscJdateResult = connection.execute(WfdiscJdateQuery)
    for WfdiscJdate in WfdiscJdateResult:
        WfdiscJdateCount = WfdiscJdate[0]
    # Integrity check on nsamp
    WfdiscNsampQuery = text(f"SELECT count(*) FROM {table_names['wfdisc']} WHERE nsamp <= 0")
    WfdiscNsampResult = connection.execute(WfdiscNsampQuery)
    for WfdiscNsamp in WfdiscNsampResult:
        WfdiscNsampCount = WfdiscNsamp[0]
    # Integrity check on samprate
    WfdiscSamprateQuery = text(f"SELECT count(*) FROM {table_names['wfdisc']} WHERE samprate <= 0")
    WfdiscSamprateResult = connection.execute(WfdiscSamprateQuery)
    for WfdiscSamprate in WfdiscSamprateResult:
        WfdiscSamprateCount = WfdiscSamprate[0]
    # Integrity check on segtype
    if CSSType == "CSS":
        WfdiscSegTypeQuery = text(
            f"SELECT count(*) FROM {table_names['wfdisc']} WHERE segtype NOT IN ('o', 'v', 's', 'd', '-')"
        )
        WfdiscSegTypeResult = connection.execute(WfdiscSegTypeQuery)
        for WfdiscSegType in WfdiscSegTypeResult:
            WfdiscSegTypeCount = WfdiscSegType[0]
    else:
        WfdiscSegTypeQuery = text(
            f"SELECT count(*) FROM {table_names['wfdisc']} WHERE segtype NOT IN ('o', 'v', 's', 'd', 'c', 'f', 'g', 'A', 'V', 'D', 'n', 't', 'u', 'x', '-')"
        )
        WfdiscSegTypeResult = connection.execute(WfdiscSegTypeQuery)
        for WfdiscSegType in WfdiscSegTypeResult:
            WfdiscSegTypeCount = WfdiscSegType[0]
    # Integrity check on time
    WfdiscTimeQuery = text(f"SELECT count(*) FROM {table_names['wfdisc']} WHERE time < -9999999999.999")
    WfdiscTimeResult = connection.execute(WfdiscTimeQuery)
    for WfdiscTime in WfdiscTimeResult:
        WfdiscTimeCount = WfdiscTime[0]
    # Integrity check on wfid
    WfdiscWfidQuery = text(f"SELECT count(*) FROM {table_names['wfdisc']} WHERE wfid <= 0")
    WfdiscWfidResult = connection.execute(WfdiscWfidQuery)
    for WfdiscWfid in WfdiscWfidResult:
        WfdiscWfidCount = WfdiscWfid[0]
    # Integrity check on datatype
    if CSSType == "CSS":
        WfdiscDatatypeQuery = text(
            f"SELECT count(*) FROM {table_names['wfdisc']} WHERE datatype NOT IN "
            "('a0', 'b0', 'c0', 'a#', 'b#', 'c#', 't4', 't8', 's4', 's2', 'f4', 'f8', 'i4', 'i2', 'g2', '-')"
        )
        WfdiscDatatypeResult = connection.execute(WfdiscDatatypeQuery)
        for WfdiscDatatype in WfdiscDatatypeResult:
            WfdiscDatatypeCount = WfdiscDatatype[0]
    else:
        WfdiscDatatypeQuery = text(
            f"SELECT count(*) FROM {table_names['wfdisc']} WHERE datatype NOT IN "
            "('a0', 'b0', 'c0', 'a#', 'b#', 'c#', 't4', 't8', 's4', 's2', 's3', 'f4', 'f8', 'i4', 'i2', 'e1', 'e#', 'g2', '-')"
        )
        WfdiscDatatypeResult = connection.execute(WfdiscDatatypeQuery)
        for WfdiscDatatype in WfdiscDatatypeResult:
            WfdiscDatatypeCount = WfdiscDatatype[0]

    # Add all results to dictionary and return
    basic_db_checks = {
        "Affiliation_EndTime_Integrity_Issues": AffEndTimeCount,
        "Affiilation_Time_Integrity_Issues": AffTimeCount,
        "Instrument_Digital_Integrity_Issues": InstDigitalCount,
        "Instrument_Inid_Integrity_Issues": InstInidCount,
        "Instrument_NCalib_Integrity_Issues": InstNcalibCount,
        "Instrument_Ncalper_Integrity_Issues": InstNcalperCount,
        "Instrument_Rsptype_Integrity_Issues": InstRsptypeCount,
        "Instrument_SampRate_Integrity_Issues": InstSampRateCount,
        "Instrument_Band_Integrity_Issues": InstBandCount,
        "Sensor_Calper_Integrity_Issues": SensorCalperCount,
        "Sensor_CalRatio_Integrity_Issues": SensorCalRatioCount,
        "Sensor_Chanid_Integrity_Issues": SensorChanidCount,
        "Sensor_EndTime_Integrity_Issues": SensorEndTimeCount,
        "Sensor_Inid_Integrity_Issues": SensorInidCount,
        "Sensor_Instant_Integrity_Issues": SensorInstantCount,
        "Sensor_Jdate_Integrity_Issues": SensorJdateCount,
        "Sensor_Time_Integrity_Issues": SensorTimeCount,
        "Sensor_Tshift_Integrity_Issues": SensorTshiftCount,
        "Site_DEast_Integrity_Issues": SiteDEastCount,
        "Site_DNorth_Integrity_Issues": SiteDNorthCount,
        "Site_Elev_Integrity_Issues": SiteElevCount,
        "Site_Lat_Integrity_Issues": SiteLatCount,
        "Site_Lon_Integrity_Issues": SiteLonCount,
        "Site_OffDate_Integrity_Issues": SiteOffDateCount,
        "Site_OnDate_Integrity_Issues": SiteOnDateCount,
        "Site_StaType_Integrity_Issues": SiteStaTypeCount,
        "Sitechan_Chanid_Integrity_Issues": SitechanChanidCount,
        "Sitechan_Ctype_Integrity_Issues": SitechanCtypeCount,
        "Sitechan_Hang_Integrity_Issues": SitechanHangCount,
        "Sitechan_Vang_Integrity_Issues": SitechanVangCount,
        "Sitechan_OffDate_Integrity_Issues": SitechanOffDateCount,
        "Sitechan_OnDate_Integrity_Issues": SitechanOndateCount,
        "Sitechan_Edepth_Integrity_Issues": SitechanEdepthCount,
        "Wfdisc_Calib_Integrity_Issues": WfiscCalibCount,
        "Wfdisc_Calper_Integrity_Issues": WfdiscCalperCount,
        "Wfdisc_Chanid_Integrity_Issues": WfdiscChanidCount,
        "Wfdisc_Clip_Integrity_Issues": WfdiscClipCount,
        "Wfdisc_Commid_Integrity_Issues": WfdiscCommidCount,
        "Wfdisc_EndTime_Integrity_Issues": WfdiscEndTimeCount,
        "Wfdisc_Foff_Integrity_Issues": WfdiscFoffCount,
        "Wfdisc_Jdate_Integrity_Issues": WfdiscJdateCount,
        "Wfdisc_Nsamp_Integrity_Issues": WfdiscNsampCount,
        "Wfdisc_Samprate_Integrity_Issues": WfdiscSamprateCount,
        "Wfdisc_SegType_Integrity_Issues": WfdiscSegTypeCount,
        "Wfdisc_Time_Integrity_Issues": WfdiscTimeCount,
        "Wfdisc_Wfid_Integrity_Issues": WfdiscWfidCount,
        "Wfdisc_Datatype_Integrity_Issues": WfdiscDatatypeCount,
    }

    imet = {
            "metric_name": "dbIntegrityMetric",
            "integrity_results": basic_db_checks
        }

    if database_config:
        database = Database(**database_config)
        database.insert_metric(imet)

    #return basic_db_checks
    return imet
