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
from obspy.core.stream import Stream
from pycheron.util.getChannel import getChannel, getChannelName, filterChannel


def test_getChannel(file_assets):

    HHE_trace = getChannel(file_assets["multi_stream"], "HHE")
    assert HHE_trace.count() > 0
    assert isinstance(HHE_trace, Stream)

    BHE_trace = getChannel(file_assets["multi_stream"], "BHE")
    assert BHE_trace.count() > 0
    assert isinstance(BHE_trace, Stream)

    H_trace = getChannel(file_assets["multi_stream"], "E")
    assert H_trace.count() > 0
    assert H_trace.count() > HHE_trace.count()
    assert H_trace.count() > BHE_trace.count()
    assert isinstance(H_trace, Stream)


def test_getChannelName():
    valid_snclq = "7A.CABN..BHE"
    valid_channel_name = getChannelName(valid_snclq)

    assert len(valid_channel_name) > 0
    assert valid_channel_name == "BHE"

    invalid_snclq = "42MAU5DFT"
    with pytest.raises(IndexError):
        getChannelName(invalid_snclq)


def test_filterChannel(file_assets):
    valid_filtered_channel = filterChannel(file_assets["BGU_HHE_stream"])

    assert valid_filtered_channel.count() > 0
    assert isinstance(valid_filtered_channel, Stream)

    # TODO: there needs to be a file with the soh_chan removed
    # in order to accomplish this, you will need an *.mseed file
    # ending in one of the values from the soh_chan list
