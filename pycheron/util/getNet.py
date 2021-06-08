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

__all__ = ["getNet", "getUniqNet", "getNetworkName"]

import numpy as np


def getNet(st, net):
    """
    Gets all traces in stream that match a specified network

    :param st: ObsPy Stream object
    :type st: obspy.core.stream.Stream
    :param net: Name of network, e.g., 7A, 6E, etc.
    :type net: str

    :return: Returns new stream object containing specific trace(s) with specified network
    :rtype: obspy.core.stream.Stream
    """

    # Grab out specified traces for station within stream object
    nets = st.select(network=net)
    return nets


def getNetworkName(snclq):
    """
    Extracts the network name from a SNCLQ id

    :param snclq: snclq id
    :type snclq: str

    :return: The name of the station
    :rtype: str
    """
    # Get network name from snclq
    netName = snclq.split(".")[0]
    return netName


def getUniqNet(st):
    """
    Obtain unique networks within a Stream object

    :param st: ObsPy Stream object
    :type st: obspy.core.stream.Stream

    :return: Return a list of unique networks within the stream object
    :rtype: list
    """

    # Create empty list variable for appending networks to
    nets = []

    # Cycle through stream object pulling out network information for each trace
    for trace in st:
        nets.append(trace.stats.network)

    # Create unique list of networks in stream
    nets = np.unique(nets)

    return nets
