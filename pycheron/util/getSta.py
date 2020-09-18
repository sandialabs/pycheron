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

__all__ = ["getSta", "getStationName", "getUniqSta"]

import numpy as np


def getSta(st, sta):
    """
    Obtain all traces in a stream that match a specified station

    :param st: Obspy stream object
    :type st: obspy.core.stream.Stream
    :param sta: Name of station, e.g., SP06, CABN
    :type sta: str

    :return: Returns a new stream object containing the specific trace(s) of the specified station
    :rtype: obspy.core.stream.Stream
    """
    # Grab out specified traces for station within stream object
    stas = st.select(station=sta)
    return stas


def getStationName(snclq):
    """
    Extracts the station name from a SNCLQ id
    :param snclq: snclq id
    :type snclq: str

    :return: The name of the station
    :rtype: str
    """
    staName = snclq.split(".")[1]
    return staName


def getUniqSta(st):
    """
    Obtain unique station information within a stream object

    :param st: Obspy stream object
    :type st: obspy.core.stream.Stream

    :return: Returns a list of unique stations from the input stream object
    :rtype: list
    """

    # Create empty list variable for appending stations to
    stas = []

    # Cycle through stream object pulling out station information for each trace
    for trace in st:
        stas.append(trace.stats.station)

    # Create unique list of stations in stream
    stas = np.unique(stas)

    return stas
