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

import obspy
from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime


def get_stream_from_file(file_name, format="MSEED"):
    return obspy.read(file_name, format=format)


def percent_diff(expected, evaluated, diff):
    return (expected == 0.0 and evaluated == 0.0) or (
        abs(expected - evaluated) / expected
    ) < diff


def compare_dicts(expected_dict, calculated_dict):
    for k, v in list(expected_dict.items()):
        if isinstance(v, (int, float)):
            assert percent_diff(v, calculated_dict[k], 0.0001)
        else:
            assert v == calculated_dict[k]


def get_stream_from_client(seis_net, channel, num, station, start_time, time_offset):
    client = Client("IRIS")
    start_time = UTCDateTime(start_time)
    end_time = start_time + time_offset
    return client.get_waveforms(seis_net, channel, num, station, start_time, end_time)


def make_dict_of_lists(list_of_dicts):
    return {
        key: [d[key] for d in list_of_dicts]
        for key, v in list(list_of_dicts[0].items())
    }
