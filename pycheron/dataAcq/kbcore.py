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
from pycheron.util.logger import Logger

__all__ = ["read_kb_core"]


def read_kb_core(kb_file, site_or_sitechan, logger=None):
    """Reads in KB Core Site file and outputs it into a pandas DF

    :param kb_file: kb core site or sitechan file
    :type kb_file: str

    :param site_or_sitechan: type of kbcore file being used
    :type site_or_sitechan: str

    :param logger: If using a logger, (you must create one using the util.logger class)
    :type logger: pycheron.util.logger.Logger

    :return: Pandas DataFrame
    :rtype: `pandas.DataFrame`

    """
    columns, _format_kb_core = _prepare_for_site_or_sitechan(site_or_sitechan)

    # kbcore site column names
    df = pd.DataFrame(columns=columns)
    with open(kb_file) as kb:
        for line in kb:
            line_as_dict = _format_kb_core(line)
            if len(line_as_dict["sta"][0]) > 0:
                df = pd.concat([df, line_as_dict], ignore_index=True)
    return df


def _prepare_for_site_or_sitechan(site_or_sitechan, logger=None):
    """Prepares values to be ingested inside read_kb_core depending on file type

    :param site_or_sitechan: type of kbcore file being used
    :type site_or_sitechan: str

    :param logger: If using a logger, (you must create one using the util.logger class)
    :type logger: pycheron.util.logger.Logger

    :return: list of column names for type of file, respective formatter function
    :rtype: `list`, function
    """
    if logger is None:
        logger = Logger(None)

    if site_or_sitechan == "site":
        columns = [
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
        _format_kb_core = _format_kb_core_site
    elif site_or_sitechan == "sitechan":
        columns = [
            "sta",
            "chan",
            "ondate",
            "chanid",
            "offdate",
            "ctype",
            "edepth",
            "hang",
            "vang",
            "descrip",
            "lddate",
        ]
        _format_kb_core = _format_kb_core_sitechan
    else:
        err_msg = "read_kb_core_site function only accepts 'site' or 'sitechan' for inputs to site_or_sitechan"
        logger.error(err_msg)
        raise ValueError(err_msg)
    return columns, _format_kb_core


def _format_kb_core_site(line):
    """Formats line of read in kb core file and separates it into columns

    :param line: kb core site file line

    :return: Pandas DataFrame
    :rtype: `pandas.DataFrame`

    """
    data = {
        "sta": line[0:7].strip(),
        "ondate": line[8:16].strip(),
        "offdate": line[17:25].strip(),
        "lat": line[26:37].strip(),
        "lon": line[38:49].strip(),
        "elev": line[50:59].strip(),
        "staname": line[59:110].strip(),
        "statype": line[110:114].strip(),
        "refsta": line[116:121].strip(),
        "dnorth": line[123:131].strip(),
        "deast": line[133:141].strip(),
        "lddate": line[143:159].strip(),
    }
    return pd.DataFrame(data, index=[0])


def _format_kb_core_sitechan(line):
    data = {
        "sta": line[0:6].strip(),
        "chan": line[7:15].strip(),
        "ondate": line[16:24].strip(),
        "chanid": line[25:33].strip(),
        "offdate": line[34:42].strip(),
        "ctype": line[43:47].strip(),
        "edepth": line[48:56].strip(),
        "hang": line[58:64].strip(),
        "vang": line[65:69].strip(),
        "descrip": line[69:120].strip(),
        "lddate": line[120:140].strip(),
    }
    return pd.DataFrame(data, index=[0])
