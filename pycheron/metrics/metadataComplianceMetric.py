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

__all__ = [
    "seedChanSpsCompliance",
    "chanOrientationCompliance",
    "seedChanConventionCheck",
    "retrieve_chan_orientation",
    "verticalChanOrientationCompliance",
    "horzChanOrientationCompliance",
    "sampleRateRespVerification",
]

from math import isclose, ceil
import numpy as np
from obspy.clients.fdsn.client import Client
from obspy.clients.iris.client import Client as ClientIRIS
from obspy import UTCDateTime
from ispaq.evalresp import evalresp as read_Evalresp
from pycheron.util.logger import Logger
from pycheron.util.getConfigFile import get_configfile
from pycheron.dataAcq.css import platform_file_path
from pycheron.dataAcq.fap import fap_reader
from obspy.io.xseed.core import _read_resp
from obspy.core.inventory.inventory import read_inventory


def seedChanSpsCompliance(st, sps_tolerance=0.01, database=None):
    """
    Checks stream object for compliance with SEED channel naming convention and sample rate. Additionally tests sample
    rate against calculated sample rate from data

    :param st: obspy stream object
    :type st: obspy.core.stream.Stream
    :param sps_tolerance: float value for allowed tolerance between metadata sample rate and actual sample rate
                          calculated from data. Also used for allowable tolerance when comparing to band code sps
                          ranges (default = 0.01 or 1%)
    :type sps_tolerance: float
    :param database: Pycheron database to ingest metrics into
    :type database: pycheron.db.sqllite_db.Database

    :return:

    list of dictionaries for each trace in a stream with the following keys and types

                * is_chan_sps_Seedcompliant (`bool`)
                * does_sps_match_data (`bool`)
                * snclq (`str`)
                * start_time (`str`)
                * end_time (`str`)
                * metric_name (`str`)

    where:
    is_chan_sps_Seedcompliant: specifies whether the bandcode (i.e., first letter of the channel name) matches
                               the sps defined in the SEED Channel Naming convention specified in Appendix A of the
                               SEED Manual. 0 means it does, 1 means it doesn't

    does_sps_match_data: specifies whether the calculated sample rate, based on calculate npts and the assumption of
                         "relatively" constant sps in the trace, matches the sample rate specified in the metadata

    snclq: station network channel location quality metadata information from the trace

    start_time: start time of the trace

    end_time: end time of the trace

    metric_name: metric name

    :rtype: list of dictionaries

    """

    # Initialize output list and masks. Masks are initialized to None in the case where the user turns off masking
    SeedChanComp = []

    # Loop through traces in stream object and pull out channel and sample rate
    for tr in st:
        chan = tr.stats.channel
        bandcode = chan[0].upper()
        sps = tr.stats.sampling_rate
        sampRatePair = seedChanConventionCheck(bandcode)

        # First check that the sample rate in the metadata matches the sample rate in the actual data.
        # Let's determine the number of points in the time series data. We won't assume that the npts in the trace
        # object is necessarily correct either. It's rare I have found this to not be the case, but not hard to
        # calculate this ourselves either.
        npts = len(tr.times())
        # Assume sample rate is constant based on the fact that Obspy enforces that
        calcSps = round(npts / (tr.stats.endtime - tr.stats.starttime))

        # If value is 1% of the sps in the metadata then there isn't an issue
        if np.abs(calcSps - sps) < sps_tolerance:
            dataSps = 1
        else:
            dataSps = 0

        # If actual sample rate matches the metadata sample rate, do the following
        if dataSps == 1:
            # Check sampRatePair against sps and bandcode to see if they match.
            # If they do, channel is in compliance with Seed Band Code naming convention (i.e., boolean value = 0),
            # otherwise it is not (i.e., boolean value = 1)
            if bandcode == "M":
                if sps > sampRatePair[0] and sps < sampRatePair[1]:
                    chanBandCode = 1
                else:
                    chanBandCode = 0

            elif bandcode == "L":
                if isclose(sps, sampRatePair, abs_tol=sps_tolerance):
                    chanBandCode = 1
                else:
                    chanBandCode = 0

            elif bandcode == "V":
                if isclose(sps, sampRatePair, abs_tol=sps_tolerance):
                    chanBandCode = 1
                else:
                    chanBandCode = 0
            elif bandcode == "U":
                if isclose(sps, sampRatePair, abs_tol=sps_tolerance):
                    chanBandCode = 1
                else:
                    chanBandCode = 0

            elif bandcode == "Q":
                if sps < sampRatePair:
                    chanBandCode = 1
                else:
                    chanBandCode = 0
            # Should probably think about a better thing to do there for 'A' and 'O' bandcodes
            elif bandcode == "A":
                chanBandCode = "N/A"
            elif bandcode == "O":
                chanBandCode = "N/A"
            else:
                if sps >= sampRatePair[0] and sps < sampRatePair[1]:
                    chanBandCode = 1
                else:
                    chanBandCode = 0
        # This means that the sample rate does not match the metadata, but it may still match the actual Band Code.
        # This block uses the calculated sample rate and compares it to the metadata bandcode to see if there still is
        # a match. May want to just have this be chanBandCode = 0.
        else:
            # Check sampRatePair against sps and bandcode to see if they match.
            # If they do, channel is in compliance with Seed Band Code naming convention (i.e., boolean value = 0),
            # otherwise it is not (i.e., boolean value = 1)
            if bandcode == "M":
                if calcSps > sampRatePair[0] and calcSps < sampRatePair[1]:
                    chanBandCode = 1
                else:
                    chanBandCode = 0

            elif bandcode == "L":
                if isclose(calcSps, sampRatePair, abs_tol=sps_tolerance):
                    chanBandCode = 1
                else:
                    chanBandCode = 0

            elif bandcode == "V":
                if isclose(calcSps, sampRatePair, abs_tol=sps_tolerance):
                    chanBandCode = 1
                else:
                    chanBandCode = 0
            elif bandcode == "U":
                if isclose(calcSps, sampRatePair, abs_tol=sps_tolerance):
                    chanBandCode = 1
                else:
                    chanBandCode = 0

            elif bandcode == "Q":
                if calcSps < sampRatePair:
                    chanBandCode = 1
                else:
                    chanBandCode = 0
            # Should probably think about a better thing to do there for 'A' and 'O' bandcodes
            elif bandcode == "A":
                chanBandCode = "N/A"
            elif bandcode == "O":
                chanBandCode = "N/A"
            else:
                if calcSps >= sampRatePair[0] and calcSps < sampRatePair[1]:
                    chanBandCode = 1
                else:
                    chanBandCode = 0

        data = {
            "is_chan_sps_Seedcompliant": chanBandCode,
            "does_sps_match_data": dataSps,
            "snclq": tr.get_id(),
            "start_time": tr.stats.starttime.isoformat(),
            "end_time": tr.stats.endtime.isoformat(),
            "metric_name": "metadataComplianceMetric",
            "metric_subname": "seedChanSpsCompliance",
        }
        SeedChanComp.append(data)
    
    if database is not None:
        database.insert_metric(SeedChanComp)
    
    return SeedChanComp


def chanOrientationCompliance(st, inv=None, angle_tolerance=0.01, iris_compatible=True, database=None):
    """

    For each trace in a stream object, check that the orientation (azimuth & dip) agree with the assumed channel
    orientation based on the channel name

    :param st: obspy stream object
    :type st: obspy.core.stream.Stream
    :param inv: obspy inventory object
    :type inv: obspy.core inventory object
    :param angle_tolerance: float value for allowed tolerance between metadata azimuth/dip and prescribed azimuth/dip
                          given channel name (default = 0.01 or 1%).  Also used to for orthogonality tolerance.
    :type angle_tolerance: float
    :param iris_compatible: flag used to let function know if data can be fetched from IRIS
    :type iris_compatible: boolean
    :param database: Pycheron database to ingest metrics into
    :type database: pycheron.db.sqllite_db.Database

    :return: list of dictionaries containing the following items for each obspy trace:
             * is_chan_orientation_compliant (boolean): 1 if both dip and azimuth are compliant;
                                                        0 if one or both are not
             * snclq: sta, net, chan, loc, quality of trace
             * start_time: start time of trace
             * end_time: end time of trace
             * metric_name: metric name, in this case ChanOrientationCompliance
    :rtype: list of dictionaries

    """
    # Cycle through traces in stream object and check for orientation compliance
    # Create list for horizontal components in the stream to iterate through after
    horzChan = []
    ChanOrientComp = []
    for tr in st:
        if inv is None:
            inv = _iris_compat_check(tr, iris_compatible)

        # First retrieve the azimuth and dip for the trace
        azimuth, dip = retrieve_chan_orientation(tr, inv)

        # Next check the vertical 'Z' component channel orientation.
        # We want to ensure that the abs(dip) is 90 and azimuth is 0.
        if azimuth is None or dip is None:
            # Should this evaluate to false or something else? Technically we don't know, so maybe the safe assumption
            # is it doesn't, or error on the side of caution. This could also be the case where the inventory object
            # provided doesn't match the stream object data
            is_chan_orientation_compliant = 0
        if tr.stats.component == "Z":
            if isclose(abs(azimuth), 0.0, abs_tol=angle_tolerance) and isclose(abs(dip), 90.0, abs_tol=angle_tolerance):
                is_chan_orientation_compliant = 1
            else:
                is_chan_orientation_compliant = 0
            data = {
                "is_chan_orientation_compliant": is_chan_orientation_compliant,
                "snclq": tr.get_id(),
                "start_time": tr.stats.starttime.isoformat(),
                "end_time": tr.stats.endtime.isoformat(),
                "metric_name": "metadataComplianceMetric",
                "metric_subname": "ChanOrientationCompliance",
            }
            ChanOrientComp.append(data)
        else:
            horzChan.append(
                (
                    tr.id,
                    azimuth,
                    dip,
                    tr.stats.starttime.isoformat(),
                    tr.stats.endtime.isoformat(),
                )
            )

    # Iterate through the list of horizontal channels comparing pairs that have matching ids and then determine whether
    # their dip is 0 and they are orthogonal to one another
    # First check that the network, sta, loc, chan[0:2] are the same, because these are the ones that we actually want
    # to.
    # To start, first grab only the trace ids to compare
    traHorzlist = [i[0] for i in horzChan]

    if len(traHorzlist) == 1:
        # TODO T: replace with logger
        print("chanOrientationCompliance(): is_chan_orientation_compliant did not check orthogonality because only one trace was provided by the user")
        
        # Check dip of original list elements
        if isclose(horzChan[0][2], 0.0, abs_tol=angle_tolerance):
            is_chan_dip_tr1 = 1
        else:
            is_chan_dip_tr1 = 0
        
        if horzChan[0][0][-1] == "E":
            if isclose(horzChan[0][1], 90.0, abs_tol=angle_tolerance):
                is_chan_az_tr1 = 1
            else:
                is_chan_az_tr1 = 0
        elif horzChan[0][0][-1] == "N":
            if isclose(horzChan[0][1], 0.0, abs_tol=angle_tolerance):
                is_chan_az_tr1 = 1
            else:
                is_chan_az_tr1 = 0
        else:
            is_chan_az_tr1 = 1

        # Now check everything to give the final verdict on compliance
        if is_chan_dip_tr1 == 1 and is_chan_az_tr1 == 1:
            is_chan_orientation_compliant1 = 1
        else:
            is_chan_orientation_compliant1 = 0

        # Fill in data for trace and append to larger list
        data = {
            "is_chan_orientation_compliant": is_chan_orientation_compliant1,
            "snclq": horzChan[0][0],
            "start_time": horzChan[0][3],
            "end_time": horzChan[0][4],
            "metric_name": "metadataComplianceMetric",
            "metric_subname": "ChanOrientationCompliance",
        }
        ChanOrientComp.append(data)

    elif len(traHorzlist) > 1:
        # Set these for channels that are '1' and '2'
        # Now iterate through list and check dip, orthogonality
        for i, (curr, next) in enumerate(zip(traHorzlist, traHorzlist[1:])):
            # If net, sta, loc, chan[0:2] are the same, then verify the dip of each (which we could've really done before)
            # and then that they are orthogonal
            if curr[:-1] == next[:-1]:
                # Check dip of original list elements
                # Trace 1
                if isclose(horzChan[i][2], 0.0, abs_tol=angle_tolerance):
                    is_chan_dip_tr1 = 1
                else:
                    is_chan_dip_tr1 = 0
                # Trace 2
                if isclose(horzChan[i + 1][2], 0.0, abs_tol=angle_tolerance):
                    is_chan_dip_tr2 = 1
                else:
                    is_chan_dip_tr2 = 0

                # If E, N component, check azimuth. Technically, if these are the expected azimuth, its proof of
                # orthogonality but we will check orthogonality anyway, as might has '1' and '2' channels as well
                # Trace 1. Technically could probably remove much of this section and just check the orthogonality directly,
                # though good for the E, N to actually verify its the standard definition
                if horzChan[i][0][-1] == "E":
                    if isclose(horzChan[i][1], 90.0, abs_tol=angle_tolerance):
                        is_chan_az_tr1 = 1
                    else:
                        is_chan_az_tr1 = 0
                elif horzChan[i][0][-1] == "N":
                    if isclose(horzChan[i][1], 0.0, abs_tol=angle_tolerance):
                        is_chan_az_tr1 = 1
                    else:
                        is_chan_az_tr1 = 0
                # We can't check azimuth here since it can be anything for '1' and '2' channels
                else:
                    is_chan_az_tr1 = 1

                # Trace 2
                if horzChan[i + 1][0][-1] == "E":
                    if isclose(horzChan[i + 1][1], 90.0, abs_tol=angle_tolerance):
                        is_chan_az_tr2 = 1
                    else:
                        is_chan_az_tr2 = 0
                elif horzChan[i + 1][0][-1] == "N":
                    if isclose(horzChan[i + 1][1], 0.0, abs_tol=angle_tolerance):
                        is_chan_az_tr2 = 1
                    else:
                        is_chan_az_tr2 = 0
                # We can't check azimuth here since it can be anything for '1' and '2' channels
                else:
                    is_chan_az_tr2 = 1

                # Now check orthogonality. Assume azimuth is right-handed coordinate system
                if abs(abs(horzChan[i][1] - horzChan[i + 1][1] - 90.0)) < angle_tolerance:
                    is_chan_orientation_orthog_compliant = 1
                else:
                    is_chan_orientation_orthog_compliant = 0

                # Now check everything to give the final verdict on compliance
                if is_chan_dip_tr1 == 1 and is_chan_az_tr1 == 1 and is_chan_orientation_orthog_compliant == 1:
                    is_chan_orientation_compliant1 = 1
                else:
                    is_chan_orientation_compliant1 = 0
                if is_chan_dip_tr2 == 1 and is_chan_az_tr2 == 1 and is_chan_orientation_orthog_compliant == 1:
                    is_chan_orientation_compliant2 = 1
                else:
                    is_chan_orientation_compliant2 = 0

                # Fill in data for both traces and append to larger list
                # Trace 1
                data = {
                    "is_chan_orientation_compliant": is_chan_orientation_compliant1,
                    "snclq": horzChan[i][0],
                    "start_time": horzChan[i][3],
                    "end_time": horzChan[i][4],
                    "metric_name": "metadataComplianceMetric",
                    "metric_subname": "ChanOrientationCompliance",
                }
                ChanOrientComp.append(data)
                # Trace 2
                data = {
                    "is_chan_orientation_compliant": is_chan_orientation_compliant2,
                    "snclq": horzChan[i + 1][0],
                    "start_time": horzChan[i + 1][3],
                    "end_time": horzChan[i + 1][4],
                    "metric_name": "metadataComplianceMetric",
                    "metric_subname": "ChanOrientationCompliance",
                }
                ChanOrientComp.append(data)
    
    if database is not None:
        database.insert_metric(ChanOrientComp)

    return ChanOrientComp


def seedChanConventionCheck(bandcode):
    """
    Function to return the provided band code sample rate pair based on the SEED band code provided
    :param bandcode: SEED bandcode from a tr.stats.channel, e.g., 'B' for broadband
    :type bandcode: string

    :return: returns the sample rate pair (or single number) based on the SEED band code provided
    :rtype: either a tuple if it's a pair or an integer
    """

    # Given the channel band code return the sampling rate pair
    # Create a dictionary of tuple values (or single values for L, V, U) to compare to the band code from the channel
    # name. Dictionary values are based on SEED Manual Appendix A: Channel Naming Convention

    bandcode_dict = {
        "F": (1000, 5000),
        "G": (1000, 5000),
        "D": (250, 1000),
        "C": (250, 1000),
        "E": (80, 250),
        "S": (10, 80),
        "H": (80, 250),
        "B": (10, 80),
        "M": (1, 10),
        "L": 1,
        "V": 0.1,
        "U": 0.01,
        "R": (0.0001, 0.001),
        "P": (0.00001, 0.0001),
        "T": (0.000001, 0.00001),
        "Q": 0.000001,
    }

    # Iterate through the keys to find the code that matches the bandcode
    sampRatePair = [val for key, val in bandcode_dict.items() if bandcode in key][0]

    return sampRatePair


def retrieve_chan_orientation(tr, inv):
    """
    Obtain the dip and azimuth of the specified trace from the provided inventory object

    :param tr: obspy trace object
    :type tr: obspy.core trace object
    :param inv: obspy inventory object
    :type inv: obspy.core inventory object

    :return: azimuth and dip of trace
    :rtype: obspy.core.inventory.util.Azimuth and obspy.core.inventory.util.Dip
    """
    # Need to ensure the trace is in the inventory object, otherwise error and set azimuth and dip to None
    try:
        # Get orientation for the trace given the provided start time
        orientn = inv.get_orientation(tr.id, tr.stats.starttime)
        azimuth = orientn["azimuth"]
        dip = orientn["dip"]
    except Exception as errmsg:
        print("Unable to retrieve azimuth and dip for %s due to: %s" % (tr.id, errmsg))
        azimuth = None
        dip = None

    return azimuth, dip


def verticalChanOrientationCompliance(tr, inv=None, angle_tolerance=0.01, iris_compatible=True, database=None):
    """

    For vertical channels only, check that the orientation of the azimuth (0 deg) and dip (-90 deg) are compliant with
    SEED convention. The individual trace and inv option could be beneficial for users that have their own response
    files

    :param tr: obspy trace object
    :type tr: obspy.core.trace object
    :param inv: obspy inventory object
    :type inv: obspy.core inventory object
    :param angle_tolerance: float value for allowed tolerance between metadata azimuth/dip and prescribed azimuth/dip
                          given channel name (default = 0.01 or 1%).
    :type angle_tolerance: float
    :param iris_compatible: flag used to let function know if data can be fetched from IRIS
    :type iris_compatible: boolean
    :param database: Pycheron database to ingest metrics into
    :type database: pycheron.db.sqllite_db.Database

    :return: dictonary containing the following items for each obspy trace:
             * is_vert_chan_orientation_compliant (boolean): 1 if both dip and azimuth are compliant;
                                                             0 if one or both are not
             * snclq: sta, net, chan, loc, quality of trace
             * start_time: start time of trace
             * end_time: end time of trace
             * metric_name: metric name, in this case verticalChanOrientationCompliance
    :rtype: dictionary

    """
    # First verify that inv is not None, if it is, throw an error
    if inv is None:
        inv = _iris_compat_check(tr, iris_compatible)
    # Ensure user only provided Z channel
    if tr.stats.component == "Z":
        # First retrieve the azimuth and dip for the vertical trace
        azimuth, dip = retrieve_chan_orientation(tr, inv)
        # Then check compliance
        if isclose(abs(azimuth), 0.0, abs_tol=angle_tolerance) and isclose(abs(dip), 90.0, abs_tol=angle_tolerance):
            is_vert_chan_orientation_compliant = 1
        else:
            is_vert_chan_orientation_compliant = 0
        data = {
            "is_vert_chan_orientation_compliant": is_vert_chan_orientation_compliant,
            "snclq": tr.get_id(),
            "start_time": tr.stats.starttime.isoformat(),
            "end_time": tr.stats.endtime.isoformat(),
            "metric_name": "metadataComplianceMetric",
            "metric_subname": "verticalChanOrientationCompliance",
        }
    else:
        raise ValueError(
            "The provided channel was %s, which is not a vertical (Z) component. "
            "verticalChanOrientationCompliance checks azimuth and dip of vertical channels only." % tr.id
        )
    
    if database is not None:
        database.insert_metric([data])

    return [data]


def horzChanOrientationCompliance(tr1, tr2, inv1=None, inv2=None, angle_tolerance=0.01, iris_compatible=True, database=None):
    """

    For horizontal channels only, check that the orientation of the azimuth and dip are compliant with
    SEED convention. Inv1 and Inv2 may be the same but we don't assume they have to be. This could be useful for cases
    where user is provided their own response files

    :param tr1: obspy trace object
    :type tr1: obspy.core.trace object
    :param inv: obspy inventory object
    :type inv: obspy.core inventory object
    :param tr2: obspy trace object
    :type tr2: obspy.core.trace object
    :param inv2: obspy inventory object
    :type inv2: obspy.core inventory object
    :param angle_tolerance: float value for allowed tolerance between metadata azimuth/dip and prescribed azimuth/dip
                          given channel name (default = 0.01 or 1%). Also used to for orthogonality tolerance.
    :type angle_tolerance: float
    :param iris_compatible: flag used to let function know if data can be fetched from IRIS
    :type iris_compatible: boolean
    :param database: Pycheron database to ingest metrics into
    :type database: pycheron.db.sqllite_db.Database
    
    :return: dictionary containing the following items for each obspy trace:
             * is_horz_chan_orientation_tr1 (boolean): 1 if both dip and azimuth are compliant; 0 if one or both are not
             * snclq_tr1: sta, net, chan, loc, quality of trace 1
             * start_time_tr1: start time of trace 1
             * end_time_tr1: end time of trace 1
             * is_horz_chan_orientation_tr2 (boolean): 1 if both dip and azimuth are compliant; 0 if one or both are not
             * snclq_tr2: sta, net, chan, loc, quality of trace 2
             * start_time_tr2: start time of trace 2
             * end_time_tr2: end time of trace 2
             * metric_name: metric name, in this case horzChanOrientationCompliance
    :rtype: dictionary

    """
    # First verify that inv1 and inv2 are not None
    if inv1 is None:
        inv1 = _iris_compat_check(tr1, iris_compatible)
    if inv2 is None:
        inv2 = _iris_compat_check(tr2, iris_compatible)

    # Ensure that network, station, loc, and first two letters of chan match, otherwise don't proceed
    if tr1.id[:-1] != tr2.id[:-1]:
        raise ValueError(
            "Trace 1 network, station, loc code, and channel [0:2] should match Trace 2, however the"
            "following traces %s, %s were provided" % (tr1.id, tr2.id)
        )

    # Error if user provided Z channel
    if tr1.stats.component == "Z" or tr2.stats.component == "Z":
        raise ValueError(
            "One or both of the provided channels %s, %s were a vertical component, which is not a "
            "horizontal component. horzChanOrientationCompliance checks azimuth and dip of horizontal "
            "channels only." % (tr1.id, tr2.id)
        )

    # First retrieve the azimuth and dip for both traces
    azimuth1, dip1 = retrieve_chan_orientation(tr1, inv1)
    azimuth2, dip2 = retrieve_chan_orientation(tr2, inv2)

    # Check dip for both traces is equal to 0
    # Trace 1
    if isclose(dip1, 0.0, abs_tol=angle_tolerance):
        is_chan_dip_tr1 = 1
    else:
        is_chan_dip_tr1 = 0
    # Trace 2
    if isclose(dip2, 0.0, abs_tol=angle_tolerance):
        is_chan_dip_tr2 = 1
    else:
        is_chan_dip_tr2 = 0

    # If E, N component, check azimuth. Technically, if these are the expected azimuth, its proof of
    # orthogonality but we will check orthogonality anyway, as might has '1' and '2' channels as well
    # Trace 1. Technically could probably remove much of this section and just check the orthogonality directly,
    # though good for the E, N to actually verify its the standard definition
    # Trace 1
    if tr1.stats.component == "E":
        if isclose(azimuth1, 90.0, abs_tol=angle_tolerance):
            is_chan_az_tr1 = 1
        else:
            is_chan_az_tr1 = 0
    elif tr1.stats.component == "N":
        if isclose(azimuth1, 0.0, abs_tol=angle_tolerance):
            is_chan_az_tr1 = 1
        else:
            is_chan_az_tr1 = 0
    # We can't check azimuth here since it can be anything for '1' and '2' channels
    else:
        is_chan_az_tr1 = 1

    # Trace 2
    if tr2.stats.component == "E":
        if isclose(azimuth2, 90.0, abs_tol=angle_tolerance):
            is_chan_az_tr2 = 1
        else:
            is_chan_az_tr2 = 0
    elif tr2.stats.component == "N":
        if isclose(azimuth2, 0.0, abs_tol=angle_tolerance):
            is_chan_az_tr2 = 1
        else:
            is_chan_az_tr2 = 0
    # We can't check azimuth here since it can be anything for '1' and '2' channels
    else:
        is_chan_az_tr2 = 1

    # Now check orthogonality. Assume azimuth is right-handed coordinate system
    if abs(abs(azimuth1 - azimuth2 - 90.0)) < angle_tolerance:
        is_chan_orientation_orthog_compliant = 1
    else:
        is_chan_orientation_orthog_compliant = 0

    # Now provide the final verdict of compliance for both traces
    if is_chan_dip_tr1 == 1 and is_chan_az_tr1 == 1 and is_chan_orientation_orthog_compliant == 1:
        is_horz_chan_orientation_compliant1 = 1
    else:
        is_horz_chan_orientation_compliant1 = 0
    if is_chan_dip_tr2 == 1 and is_chan_az_tr2 == 1 and is_chan_orientation_orthog_compliant == 1:
        is_horz_chan_orientation_compliant2 = 1
    else:
        is_horz_chan_orientation_compliant2 = 0

    # Output data to dictionary
    data = {
        "is_horz_chan_orientation_compliant_tr1": is_horz_chan_orientation_compliant1,
        "snclq_tr1": tr1.get_id(),
        "start_time_tr1": tr1.stats.starttime.isoformat(),
        "end_time_tr1": tr1.stats.endtime.isoformat(),
        "is_horz_chan_orientation_compliant_tr2": is_horz_chan_orientation_compliant2,
        "snclq_tr2": tr2.get_id(),
        "start_time_tr2": tr2.stats.starttime.isoformat(),
        "end_time_tr2": tr2.stats.endtime.isoformat(),
        "metric_name": "metadataComplianceMetric",
        "metric_subname": "horzChanOrientationCompliance",
    }

    if database is not None:
        database.insert_metric([data])

    return [data]


def sampleRateRespVerification(st, perc_tol=0.15, norm_freq=None, iris_compatible=True, database=None, logger=None):
    """
    Sample Rate consistency check that compares the sample rate from a stream object to the sample rate obtained from
    the high frequency corner of the trace's (i.e., channel's) amplitude response
    (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf).
    Returns a boolean value, where 0 indicates that the miniSEED sample rate and the response derived sample rate
    agree within perc_tol, or 1 if they disagree.

    :param st: ObsPy Stream object
    :type st: obspy.core.stream.Stream
    :param perc_tol: maximum % difference allowed for stream and response derived sample rates to be considered
                     "close" in value (default = 0.15, or 15%)
    :type: float
    :param norm_freq: normalization frequency that the instrument sensitivity or dataless stage 0 sensitivity is valid
                      (default = None)
    :param iris_compatible: flag used to let function know if data can be fetched from IRIS
    :type iris_compatible: boolean
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger

    :return: Returns list of dictionaries that includes the following for each trace within the stream object:

            * snclq (`str`): station network channel location quality code information for specfic trace
            * start_time (`str`): start time of trace
            * end_time (`str`): end time of trace
            * sample_rate_resp (`boolean`): sample_rate_resp, boolean value, where 0 indicates that the miniSEED sps
              and response derived sample rate agree within perc_tol, or 1 if they disagree
            * metric_name (`str`): sampleRateRespVerification string
    :rtype: list

    * Code originally ported from IRISMustangMetrics R Cran Package
    (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, March 20). CRAN-Package IRISMustangMetrics.
    The Comprehensive R Archive Network. Retrieved from
    https://cran.r-project.org/web/packages/IRISMustangMetrics/index.html) and augmented and adapted for use within
    Pycheron

    **Algorithm steps:**
    Summarized from steps below and paraphrased from
    (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)
    #. Retrieves instrument response from IRIS or from user provided response file using frequencies that are
       one decade below the normalization frequency and one decade above the miniSEED sampling rate. Also gathers
       instrument sensitivity from inventory object to use as normalization frequency if not provided by user as input.
    #. Obtains the difference in amplitude values, normalized for frequency spacing, so that the first steep rolloff
       can be determined
    #. Scans through the amplitude differences to find the first steep rolloff, maximum difference in the rolloff
       (arg.min)
    #. The corresponding frequency value at the first steep rolloff is utilized as the corner frequency and multiplied
       by two to provide the empirical response derived sample rate
    #. The empirical response derived sample rate is then compared to the miniSEED sample rate, and if it is within the
       perc_tol defined by the user (default is 15 (same as IRIS) as there can be variations across instruments) then
       sample_rate_resp = 0 indicating the sample rates agree. Otherwise sample_rate_resp = 1 if they disagree

    """
    # Set up IRIS Client
    IRIS = ClientIRIS(timeout=30)
    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Set up output list. Will be list of dictionaries
    d = []

    # Set up appropriate metadata for requests, if needed
    for tr in st:
        starttime = tr.stats.starttime
        endtime = tr.stats.endtime
        mseed_sps = tr.stats.sampling_rate
        net = tr.stats.network
        sta = tr.stats.station
        loc = tr.stats.location
        chan = tr.stats.channel

        if iris_compatible:
            # Check if normalization frequency is null, if so grab it from IRIS
            # TODO grab this from a user supplied file too
            client = Client("IRIS")
            inv = client.get_stations(
                network=net,
                station=sta,
                location=loc,
                starttime=starttime,
                endtime=endtime,
                level="channel",
            )
            netInv = inv[0]
            staInv = netInv[0]
            chanInv = staInv[0]
            norm_freq = chanInv.response.instrument_sensitivity.frequency

            # After we've grabbed the info above, either from the user provided input or from IRIS
            # ensure normalization frequency is not zero and is not null
            if norm_freq is None or norm_freq == 0:
                logger.error(f"metadataComplianceMetric(): sampleRateRespVerification error, normalization frequency is None or 0: norm_freq={norm_freq}")
                return

            minFreq, maxFreq, nFreq = _calc_maxf_maxf_nf(norm_freq, mseed_sps)
            units = "def"

            # Check if min frequency greater than or equal to maximum frequency
            if minFreq >= maxFreq:
                logger.error(
                    "metadataComplianceMetric(): sampleRateRespVerification error, calculated minimum frequency "
                    "(from normalization frequnecy) is greater than or equal to calculated maximum frequency "
                    "(from miniSEED sample rate)"
                )
                return

            # Obtain response information
            # Get instrument response from IRIS webservice in the case where evalresp = None
            try:
                resp = IRIS.evalresp(
                    net,
                    sta,
                    loc,
                    chan,
                    starttime,
                    minFreq,
                    maxFreq,
                    nFreq,
                    units,
                    output="fap",
                )
                amp = resp[:, 1]
                freq = resp[:, 0]
            except Exception as e:
                raise ValueError(f"Failed to retrieve IRIS evalresp data in metadataComplianceMetric: {e}")
        
        # Otherwise read instrument derived file, either FAP or EVRESP
        else:
            evalresp = get_configfile(tr, config_type="response")
            e_filtered = evalresp[chan]
            respfile = platform_file_path(filename=e_filtered[1], dirname=e_filtered[0])

            if e_filtered[2].upper() == "FAP":
                if norm_freq is None or norm_freq == 0:
                    logger.error(f"metadataComplianceMetric(): sampleRateRespVerification error, normalization frequency is None or 0: norm_freq={norm_freq}")
                    return
                resp = fap_reader(respfile)

            elif e_filtered[2].upper() == "PAZFIR":
                logger.error("metadataComplianceMetric(): Current supported file types are FAP and EVRESP")
                raise ValueError("Current supported file types are FAP and EVRESP")

            elif e_filtered[2].upper() == "EVRESP":
                inv_read = _read_resp(respfile)
                norm_freq = inv_read[0][0][0].response.instrument_sensitivity.frequency

                if norm_freq is None or norm_freq == 0:
                    logger.error(f"metadataComplianceMetric(): sampleRateRespVerification error, normalization frequency is None or 0: norm_freq={norm_freq}")
                    return

                minFreq, maxFreq, nFreq = _calc_maxf_maxf_nf(norm_freq, mseed_sps)
                units = "def"

                # Check if min frequency greater than or equal to maximum frequency
                if minFreq >= maxFreq:
                    logger.error(
                        "metadataComplianceMetric(): sampleRateRespVerification error, calculated minimum frequency "
                        "(from normalization frequnecy) is greater than or equal to calculated maximum frequency "
                        "(from miniSEED sample rate)"
                    )
                    return

                try:
                    resp = read_Evalresp(
                        sfft=minFreq,
                        efft=maxFreq,
                        nfft=nFreq,
                        filename=respfile,
                        date=UTCDateTime(starttime),
                        station=sta,
                        channel=chan,
                        network=net,
                        locid=loc,
                        units=units.upper(),
                        debug=True,
                    )
                except Exception as e:
                    raise ValueError(f"Failed to read EVRESP file in metadataComplianceMetric: {e}")

            freq, amp, phase = resp[0], resp[1], resp[2]

        # Add check for if Evalresp empty, returns all zero values, or returns all null values
        # Empty response
        if len(resp) == 0:
            logger.error("metadataComplianceMetric(): sampleRateRespVerification error, empty response")
            return
        # Response all 0's
        if np.all(resp[1] == 0):
            logger.error("metadataComplianceMetric(): sampleRateRespVerification error, amplitude response all zero")
            return
        # Response None/NaNs
        if np.all(np.isnan(resp[1])) or np.all(resp[1]) is None:
            logger.error(
                "metadataComplianceMetric(): sampleRateRespVerification error, amplitude response all None/Nans"
            )
            return

        dAdf = [float("NaN")]
        # Obtain the difference in amplitude values, normalized for frequency spacing, so that we can find the first
        # steep rolloff below
        dAmp = subtract_lists(amp[1:], amp[:-1])
        dFreq = subtract_lists(freq[1:], freq[:-1])
        diffAmp = [i / j for i, j in zip(dAmp, amp[1:])]
        absFreq = [abs(i) for i in dFreq]
        diffFreq = [i / j for i, j in zip(absFreq, freq[1:])]
        dAdf.extend([i / j for i, j in zip(diffAmp, diffFreq)])

        # Set a few variables to start
        foundRollOff = False
        cornerFreq = 0

        # Loop through dAdf to find the rolloff and the corresponding frequency that is affiliated with the maximum
        # negative slope. Use -50, as empirically it implies FIR rolloffs, but not sensor poles
        # (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf)
        # Create list to append value where i < -50
        rollOff = []
        for i in range(len(dAdf)):
            if np.isnan(dAdf[i]):
                continue
            if dAdf[i] > -50 and not foundRollOff:
                continue
            if dAdf[i] <= -50:
                foundRollOff = True
                rollOff.append((freq[i], amp[i], dAdf[i]))
            if dAdf[i] > -50 and foundRollOff:
                arrRollOff = np.array(rollOff)
                cornerFreq = arrRollOff[:, 0][np.argmin(arrRollOff[:, 2])]
                break

        # Calculate empirical response derived sample rate
        resp_sps = 2 * cornerFreq

        # Determine if resp_sps falls within +/- percent tolerance (perc_tol) of the miniSEED sps
        # Inside percent tolerance
        if isclose(resp_sps, mseed_sps, rel_tol=perc_tol):
            sample_rate_resp = 0
        # Outside percent tolerance
        else:
            sample_rate_resp = 1
        # Create dictionary output object
        spsResp = {
            "snclq": tr.get_id(),
            "start_time": tr.stats.starttime.isoformat(),
            "end_time": tr.stats.endtime.isoformat(),
            "sample_rate_resp": sample_rate_resp,
            "metric_name": "metadataComplianceMetric",
            "metric_subname": "sampleRateRespVerification",
        }
        # Append trace information to list
        d.append(spsResp)

    # If database defined, insert metric information
    if database is not None:
        database.insert_metric(d)

    return d


def subtract_lists(list1, list2):
    """
    Subtracts two lists
    :param list1: list 1 to subtract from list 2
    :type: list
    :param list2: list 2 to subtract from list 1
    :type: list
    """
    return [x1 - x2 for x1, x2 in zip(list1, list2)]


def _calc_maxf_maxf_nf(norm_freq, mseed_sps):
    """
    Calculates the min/max frequency and nFreq given the normalization frequency
    :param norm_freq: normalization frequency
    :type: float
    :param mseed_sps: trace sampling rate
    :type mseed_sps: float
    """
    # Retrieve frequencies that are 1 decade below the normalization frequency and one decade above the
    # miniseed sample rate, i.e., 100 frequency samples per decade for a given response
    minFreq = norm_freq / 10
    maxFreq = 10 * mseed_sps
    nFreq = ceil(np.log10(maxFreq / minFreq) * 100)
    return (minFreq, maxFreq, nFreq)


def _iris_compat_check(tr, iris_compatible):
    """
    Checks to see if iris_compatible flag was set. If so, inventory data is fetched from IRIS.
    Otherwise, inventory file specified in inventoryfile_config.toml is used. 
    :param tr: obspy trace object
    :type tr: obspy.core trace object
    :param iris_compatible: flag used to let function know if data can be fetched from IRIS
    :type iris_compatible: boolean
    """
    # If we can reach out to IRIS to retrieve inventory information, do so
    if iris_compatible:
        client = Client("IRIS")
        inv = client.get_stations(
                network=tr.stats.network,
                station=tr.stats.station,
                starttime=tr.stats.starttime,
                endtime=tr.stats.endtime,
                level="response",
        )
    # Otherwise, use inventoryfile_config.toml. Files MUST be of type STATIONXML
    else:
        inv_loc = get_configfile(tr, config_type="inventory")
        inv_chan = inv_loc[tr.stats.channel]
        inv_path = platform_file_path(filename=inv_chan[1], dirname=inv_chan[0])
        inv = read_inventory(inv_path, format="STATIONXML")
    return inv