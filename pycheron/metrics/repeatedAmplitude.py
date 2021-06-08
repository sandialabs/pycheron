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

import os
import subprocess

import numpy as np
import pandas as pd
from pathlib2 import Path

from pycheron.util.logger import Logger

__all__ = ["_repeatedAmplitudeMetric", "repeatedAmplitudeMetric"]

import platform

if platform.system() != "Windows":
    try:
        import pycheron.repAmpsMetric as repAmpsMetric
    except ImportError:
        fpath = os.path.dirname(__file__)
        pwd = os.getcwd()
        if pwd != str(Path(fpath).parent):
            os.chdir(str(Path(fpath).parent))
        print("------------------------------------------------------")
        print("Building Fortran Library")
        print("------------------------------------------------------")
        subprocess.call(
            [
                "f2py",
                "-c",
                "-m",
                "repAmpsMetric",
                "--quiet",
                os.path.dirname(__file__) + "/repAmps/repAmps.f90",
            ]
        )
        os.chdir(pwd)
        try:
            import pycheron.repAmpsMetric as repAmpsMetric
        except ImportError:
            import repAmpsMetric


def repeatedAmplitudeMetric(st, minRep=10, generateMasks=False, logger=None, fortran=False, database=None):
    """
    Wrapper for repeatedAmplitudeMetric so entire stream will be processed

    :param st: Stream Object
    :type st: obspy.core.stream.Stream
    :param minRep: Minimum number of repeated adjacent values in a series, > 1
    :type minRep: int
    :param generateMasks: If true, generate boolean qc mask.
    :type generateMasks: bool
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger
    :param fortran: Use Fortran libs or not. If libs will not compile or on a Windows Machine, set to False
    :type fortran: bool
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database

    :return: list of dictionaries for each trace containing the following keys and types:

        * snclq (`str`)
        * count (`int`) - Number of Repeated Amplitudes in trace
        * rep_amp (`list`) - A list containing dictionaries of the Repeated Amplitude information:

          * start_time (`str`) - starttime of the repeated amplitude
          * end_time (`str`) - endtime of the repeated amplitude
          * duration (`float`) - the length (in seconds) of the repeated amplitude
          * value (`str`) - The value in the trace that is repeated

        * masks (`numpy.ndarray` or `bool`)

          * generateMasks = `True`:

            * Boolean mask: 0 (False) or 1 (True)
            * "No Mask Created": No mask was created because there were no qc issues

          * generateMasks = `False`:

            * `False`

        * metric_name (`str`)
        * start_time (`str`)
        * end_time (`str`)

    :rtype: list

    **Example**

    .. code-block:: python

        # Initialize IRIS client
        client = Client("IRIS")

        # Define start/end time and then get stream object
        starttime = UTCDateTime("2012-12-12T00:00:00.000")
        endtime = UTCDateTime("2012-12-13T00:00:00.000")
        st = client.get_waveforms("AK","GHO","","BHN",starttime, endtime)

        tr = st[0]
        minRep = 10
        repAmp = repeatedAmplitudeMetric(tr, minRep)
        print "count: ", repAmp['count']
        >>> count:  1
        print "SNCLQ: ", repAmp['snclq']
        >>> SNCLQ:  AK.GHO..BHN
        print "Masks: ", repAmp['mask']
        >>> Masks:  False

        for i in range(len(repAmp['rep_amp'])): # rep amps
            print "duration: ", repAmp['rep_amp'][i]['duration']
            print "start_time: " , repAmp['rep_amp'][i]['start_time']
            print "end_time: ", repAmp['rep_amp'][i]['end_time']
            print "value: ", repAmp['rep_amp'][i]['value']
        >>> duration:  0.2
        >>> start_time:  2012-12-12T07:45:29.468400Z
        >>> end_time:  2012-12-12T07:45:29.648400Z
        >>> value:  1292

        print "metric_name: ",  repAmp['metric_name']
        >>> metric_name:  repeatedAmplitude
        print "start_time: " , repAmp['start_time']
        >>> start_time:  2012-12-12T00:00:00.008400Z
        print "end_time: ", repAmp['end_time']
        >>> end_time:  2012-12-12T23:59:59.988400Z


    """

    # Initialize out list
    out = []
    # Loop through each trace in the stream and determine if repeated amplitude values exist
    for i in range(len(st)):
        tr = st[i]
        repAmps = _repeatedAmplitudeMetric(tr, minRep, generateMasks, logger, fortran)
        # Append repAmps to output
        out.append(repAmps)

    # If database defined, insert metric information
    if database is not None:
        database.insert_metric(out)

    return out


def _repeatedAmplitudeMetric(tr, minRep=10, generateMasks=False, logger=None, fortran=False):
    """

    Function to find repeated amplitude values

    :param tr: trace object
    :param minRep: (int) minimum number of repeated adjacent values in a series, > 1
    :param generateMasks: (bool) generate boolean qc mask. Default = False
    :param logger: (logger object) - If using a logger, (you must create one using the util.logger class)
                   (DEFAULT = None)

    :param fortran: (boolean) - Use Fortran libs or not. If libs will not compile or on a Windows Machine, set to False
                                (DEFAULT = True)

    :return: dictionary containing:
        count - Number of Repeated Amplitudes in tr
        repAmp - A list containing dictionaries of the Repeated Amplitude information:
            start_time - starttime of the repeated amplitude
            end_time - endtime of the repeated amplitude
            duration - the length (in seconds) of the repeated amplitude
            value - The value in the trace that is repeated
        masks -
            generateMasks = True:
                - Boolean mask: 0 (False) or 1 (True)
                - "No Mask Created": No mask was created because there were no qc issues
            generateMasks = False:
                - False
    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Get sample rate and trace data
    sampRate = tr.stats.sampling_rate
    x = tr.data

    # If not using fortran
    if not fortran:
        vals = []
        index = []
        # goes thru data and checks if from data point i to i + minRep they are the same.
        for i in range(minRep, len(x) - minRep):
            if x[i] == x[i + minRep]:
                if len(np.unique(x[i : i + minRep])) == 1:
                    for j in np.arange(i, i + minRep):
                        vals.append(x[j])
                        index.append(j)
    # Otherwise
    else:

        trunc, i, v, t = repAmpsMetric.repamps(x, minRep, sampRate, 200)

        index = i[0:trunc]
        vals = v[0:trunc]

    # Get unique values and initialize values/lists
    un_values = np.unique(vals)
    count = 0
    repAmps = []

    # loop through segments that have unique values
    for j in un_values:
        # Get start time
        starttime = tr.stats.starttime
        # get the index of i in the list vals
        arr = [index[i] for i, y in enumerate(vals) if y == j]
        un_arr = np.unique(arr)
        # creates breaks between segments.
        val_break = [0, len(un_arr)]
        if np.max(un_arr) - np.min(un_arr) + 1 != len(un_arr):
            for i in range(1, len(un_arr)):
                if un_arr[i] - un_arr[i - 1] != 1:
                    val_break.append(i)
            val_break.sort()
            for k in range(len(val_break) - 1):
                in_start = val_break[k]
                in_end = val_break[k + 1] - 1
                length = (un_arr[in_end] - un_arr[in_start]) / sampRate
                start = starttime + (un_arr[in_start] / sampRate)
                end = starttime + (un_arr[in_end] / sampRate)
                count += 1
                d = {
                    "start_time": start.isoformat(),
                    "end_time": end.isoformat(),
                    "value": str(j),
                    "duration": length,
                }
                repAmps.append(d)
        else:
            length = len(un_arr) / sampRate
            start = starttime + (un_arr[0] / sampRate)
            end = starttime + (un_arr[-1] / sampRate)
            count += 1

            d = {
                "start_time": start.isoformat(),
                "end_time": end.isoformat(),
                "value": str(j),
                "duration": length,
            }
            repAmps.append(d)

    # creating qc masks
    if generateMasks:
        if repAmps:
            df = pd.DataFrame(repAmps)
            df.drop("value", 1)
            df.drop("duration", 1)
            m = df.to_dict()
        else:
            m = None
    else:
        m = None

    # Generate out dictionary
    out = {
        "snclq": tr.get_id(),
        "count": count,
        "rep_amp": repAmps,
        "mask": m,
        "metric_name": "repeatedAmplitudeMetric",
        "start_time": tr.stats.starttime.isoformat(),
        "end_time": tr.stats.endtime.isoformat(),
    }

    return out
