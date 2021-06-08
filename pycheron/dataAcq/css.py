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

from obspy.core.stream import Stream
from obspy.core.trace import Trace
from obspy.core.utcdatetime import UTCDateTime
from pisces.io.readwaveform import read_waveform
import os
import platform
import pandas as pd
import numpy as np


__all__ = [
    "read_wfdisc",
    "get_wfdisc_stations",
    "css2stream",
    "wf_dir_path_create",
    "create_wfdisc_header",
    "one_stream_per_day",
    "css2streamDB",
    "platform_file_path",
]


def read_wfdisc(file):
    """Reads in wfdisc file into a Pandas DataFrame.

    :param file: wfdisc file
    :type file: str

    :return: Pandas DataFrame
    :rtype: `pandas.DataFrame`

    """
    # wfdisc column names
    columns = [
        "sta",
        "chan",
        "time",
        "wfid",
        "chanid",
        "jdate",
        "endtime",
        "nsamp",
        "samprate",
        "calib",
        "calper",
        "instype",
        "segtype",
        "datatype",
        "clip",
        "dir",
        "dfile",
        "foff",
        "commid",
        "lddate",
    ]
    return pd.read_csv(file, delim_whitespace=True, names=columns, encoding="utf-8", index_col=False)


def wf_dir_path_create(wfdisc_dir):
    """Finds first wfdisc file, appends directory it's located in to file name

    :param wfdisc_dir: the directory containing wfdisc
    :type wfdisc_dir: str

    :return: single wfdisc file path
    :rtype: str

    """
    wf_file = [os.path.join(wfdisc_dir, f) for f in os.listdir(wfdisc_dir) if f.endswith(".wfdisc")]
    return wf_file[0]


def get_wfdisc_stations(wfdisc_dir):
    """Gets all of the stations contained in wfdisc file

    :param wfdisc_dir: the directory containing wfdisc
    :type wfdisc_dir: str

    :return: list of all stations in wfdisc
    :rtype: `list`

    """
    wf_file = wf_dir_path_create(wfdisc_dir)
    wfdisc = read_wfdisc(wf_file)
    stations = np.unique(wfdisc["sta"])
    return stations


def create_wfdisc_header(wfdisc_row, network):
    """Creates header to be used in obspy Trace (Trace.stats)

    :param wfdisc_row: single row of wfdisc dataframe
    :type wfdisc_row: pandas.core.series.Series

    :return: dictionary of header info for obspy Trace
    :rtype: dict

    """
    header = {
        "sampling_rate": wfdisc_row["samprate"],
        "calib": wfdisc_row["calib"],
        "npts": wfdisc_row["nsamp"],
        "channel": wfdisc_row["chan"],
        "station": wfdisc_row["sta"],
        "network": network,
        "starttime": UTCDateTime(wfdisc_row["time"]),
        "endtime": UTCDateTime(wfdisc_row["endtime"]),
    }
    return header


def one_stream_per_day(starttimes_list, traces):
    """Takes in a list of wfdisc starttimes and traces to return a list of streams

    :param starttimes_list: list of wfdisc starttimes
    :type starttimes_list: list[obspy.core.utcdatetime.UTCDateTime]
    :param traces: obspy traces
    :type traces: from obspy.core.trace.Trace

    :return: list of obspy streams
    :rtype: list[obspy.core.trace.Stream]

    """
    st = []
    starttimes = np.unique(starttimes_list)

    for time in starttimes:
        dayTraces = [tr for tr in traces if tr.stats.starttime.isoformat()[0:10] == time]
        st.append(Stream(dayTraces))
    return st


def css2stream(wfdisc_dir, network, station, byDay=False):
    """Takes in a wfdisc file, reads it, and takes the waveform files and converts them to a stream

    :param wfdisc_dir: Directory containing wfdisc file
    :type wfdisc_dir: str
    :param network: Network which the stations in the wfdisc belong to
    :type network: str
    :param station: The station which should be extraced from file and stream created
    :type station: str
    :param byDay: Determines whether you want 1 Stream/day for station
    :type byDay: bool

    :return: Stream for certain station
    :rtype: `obspy.core.stream.Stream`

    """

    wf_file = wf_dir_path_create(wfdisc_dir)
    wfdisc = read_wfdisc(wf_file)
    # making traces
    traces = []
    starttimes_list = []
    # loop through column names
    for i, row in wfdisc.iterrows():
        if station == row["sta"]:
            header = create_wfdisc_header(row, network)
            starttimes_list.append(UTCDateTime(row["time"]).isoformat()[0:10])

            # TODO: See how this would handle a non-relative path provided in
            wfile = os.path.join(wfdisc_dir, row["dir"], row["dfile"])

            # reading waveform
            data = read_waveform(wfile, row["datatype"], row["foff"], row["nsamp"])
            tr = Trace(data, header)
            traces.append(tr)
    # if want entire station Stream
    if byDay is False:
        return Stream(traces)
    # if want 1 Stream/day for station
    else:
        return one_stream_per_day(starttimes_list, traces)


def css2streamDB(wfdisc, network, station, byDay=False):
    """Takes in a wfdisc file, reads it, and takes the waveform files and converts them to a stream

    :param wfdisc_dir: wfdisc dataframe
    :type wfdisc_dir: str
    :param network: Network which the stations in the wfdisc belong to
    :type network: str
    :param station: The station which should be extraced from file and stream created
    :type station: str
    :param byDay: Determines whether you want 1 Stream/day for station
    :type byDay: bool

    :return: Stream for certain station
    :rtype: `obspy.core.stream.Stream`

    """
    # making traces
    traces = []
    starttimes_list = []
    # loop through column names
    for i, row in wfdisc.iterrows():
        if station == row["sta"]:
            header = create_wfdisc_header(row, network)
            starttimes_list.append(UTCDateTime(row["time"]).isoformat()[0:10])
            wfile = platform_file_path(filename=row["dfile"], dirname=row["dir"])

            # reading waveform
            data = read_waveform(wfile, row["datatype"], row["foff"], row["nsamp"])
            tr = Trace(data, header)
            traces.append(tr)
    # if want entire station Stream
    if byDay is False:
        return Stream(traces)
    # if want 1 Stream/day for station
    else:
        return one_stream_per_day(starttimes_list, traces)


def platform_file_path(filename, dirname):
    """Formats path to file for various OS's

    :param filename: name of file
    :type  filename: str
    :param dirname: name of directory
    :type  dirname: str

    :returns: formatted path
    :rtype: str

    """
    # If client is Windows and DB Path is Linux or Mac
    if platform.system() == "Windows" and len(dirname.split("/")) > 1:
        windows_dir = dirname.replace("/", "\\")
        combined_file = os.path.join(windows_dir, filename)
    # If client is Linux or Mac and DB Path is Windows
    elif platform.system() != "Windows" and len(dirname.split("\\")) > 1:
        unix_dir = dirname.replace("\\", "/")
        combined_file = os.path.join(unix_dir, filename)
    else:
        combined_file = os.path.join(dirname, filename)
    return combined_file
