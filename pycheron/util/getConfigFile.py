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

from obspy import UTCDateTime
import toml
import os
from pathlib import Path

__all__ = ["get_configfile"]


def get_configfile(
    tr,
    config_type,
):
    """
    Obtains config file from either inventoryfile_config.toml or responsefile_config.toml

    :param tr: obspy trace object
    :type tr: obspy.core.trace object
    :param config_type: type of config file to be used, accepted values of 'inventory' or 'response'
    :type config_type: str

    :return: dictionary containing file info
    :rtype: dict {network: {station: {location: {channel: file info}}}}
    """
    if config_type not in ("inventory", "response"):
        raise ValueError("Invalid config_type, can not retrieve config file info for inventory or resonse data. Use either 'inventory' or 'response'")

    # Different traces have different channels
    # This means, we'll need to run every trace in a loop for each response file
    file_path = Path(os.path.realpath(__file__))
    conf_path = file_path.parent.parent.parent
    conf_file = Path(conf_path, f"{config_type}file_config.toml")
    conf = toml.load(conf_file)

    net = tr.stats.network
    sta = tr.stats.station
    loc = tr.stats.location if tr.stats.location != "" else "--"

    try:
        chan_resps = conf[net][sta][loc]
        return chan_resps
    except Exception as e:
        raise ValueError(f"Unable to read values from {config_type}file_config.toml. Could not find network, station, location, channel in config: {e}")