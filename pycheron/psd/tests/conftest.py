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
from obspy import read
from obspy.core.trace import Trace


from pycheron.psd.psdList import psdList


@pytest.fixture(scope="session")
def file_assets():
    file_info = {
        "dir": "test_data/",
        "mseeds": [
            "test_data_2011_05_01_number_1.MSEED",
            "test_data_2011_05_01_number_2.MSEED",
            "BGU_HHE_001.mseed",
            "test_data_2010_02_27_6_to_9.MSEED",
            "LHZ_Channel_Chile_earth_quake.MSEED",
            "M_chan_II_ESK_00_MHE.mseed",
            "V_channel_wave_data.MSEED",
        ],
    }

    test_streams = [
        read(file_info["dir"] + trace_file) for trace_file in file_info["mseeds"]
    ]
    file_info["test_streams"] = test_streams
    test_traces = [stream[0] for stream in test_streams]
    file_info["test_traces"] = test_traces
    return file_info


@pytest.fixture(scope="session")
def BHZ_channel_stream_and_psd_list(file_assets):
    stream_ge_one_hour = file_assets["test_streams"][3]
    return stream_ge_one_hour, psdList(stream_ge_one_hour)
