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

import numpy as np
from obspy.core.stream import Stream
from pycheron.util.getNet import getNet, getNetworkName, getUniqNet


def test_getNet(file_assets):
    # Check that the stream filters out only 7A networks
    filtered_network = getNet(file_assets["multi_stream"], "7A")

    assert filtered_network.count() > 0
    assert isinstance(filtered_network, Stream)

    # Check that the stream is empty after filtering for 7A networks
    # on UU network
    empty_filtered_network = getNet(file_assets["BGU_HHE_stream"], "7A")
    assert empty_filtered_network.count() == 0
    assert isinstance(empty_filtered_network, Stream)


def test_getNetworkName():
    # Check that, given a valid network, the name is returned
    valid_snclq = "7A.CABN..BHE"
    valid_network_name = getNetworkName(valid_snclq)

    assert len(valid_network_name) > 0
    assert valid_network_name == "7A"

    # Check that, given an empty network, an empty string is returned
    missing_network_snclq = ".CABN..BHE"
    empty_network_name = getNetworkName(missing_network_snclq)

    assert len(empty_network_name) == 0
    assert empty_network_name == ""


def test_getUniqNet(file_assets):
    # Return networks, and ensure they are unique
    unique_network = getUniqNet(file_assets["multi_stream"])

    assert isinstance(unique_network, np.ndarray)
    assert len(file_assets["mseeds"]) >= unique_network.size
    assert len(set(unique_network)) == unique_network.size
