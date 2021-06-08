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


from pycheron.dataAcq import readCSS as css
from pycheron.util.logger import Logger

wfdisc = "/Users/jbobeck/DNE18/Synthetic_QC_fresh/"

st = css.css2stream(wfdisc, "UU", "NOQ", True)

# make jitter.dat
# pd.DataFrame(st[1][2].data[4600000:5400000]).T.to_csv("jitter.dat", sep='\t', index=False, header=False)

tr = st[1][2].data[4600000:5400000]

tr = [
    2,
    3,
    4,
    60,
    60,
    60,
    60,
    60,
    60,
    60,
    60,
    60,
    60,
    60,
    60,
    60,
    3,
    2,
    1,
    5,
    4,
    3,
    2,
    5,
    5,
    1,
    2,
    4,
    30,
    30,
    30,
    30,
    30,
    6,
    3,
    4,
    6,
    7,
    8,
    0,
    1,
    5,
    6,
    3,
    5,
    3,
    1,
    3,
    1,
    3,
    4,
    2,
    20,
    25,
    22,
    23,
    27,
    30,
    31,
    29,
    31,
    30,
    30,
    31,
]
minRep = 4
# fs = 100.0
bound = 2


def repeatedAmplitudeMetric(tr, minRep=10, bound=10, generateMasks=False, logger=None, fortran=False):
    if logger is None:
        logger = Logger(None)

    fs = tr.stats.sampling_rate
    x = tr.data
    repAmps = []
    if not fortran:
        i = 0
        while i <= (len(x) - minRep):
            # whether or not to increment i
            inc = True
            count = 1
            # case 1 + 2

            if x[i] == x[i + minRep]:
                if np.all([j == x[i] for j in x[i : i + (minRep + 1)]]):  # start and endtime is the same
                    while np.all([j == x[i] for j in x[i : i + (count)]]):
                        count += 1
                    if count >= minRep:
                        type = "Flat"
                        starttime = tr.stats.starttime + (i * fs)
                        i = i + (count - 1)
                        endtime = tr.stats.starttime + (i * fs)
                        val = x[i]
                        duration = count * fs
                        d = {
                            "start_time": starttime.isoformat(),
                            "end_time": endtime.isoformat(),
                            "value": val,
                            "duration": duration,
                            "type": type,
                        }
                        repAmps.append(d)
                        inc = False
            # case 3+4
            if x[i] == x[i + minRep] or abs(x[i] - x[i + minRep]) <= bound:
                if np.all([abs(j - x[i]) <= bound for j in x[i : i + (minRep + 1)]]):
                    while np.all([abs(j - x[i]) < bound for j in x[i : i + (count)]]) and i + count < len(x):
                        count += 1

                    if count >= minRep:
                        type = "Jitter"
                        starttime = tr.stats.starttime + (i * fs)
                        i = i + (count - 2)
                        endtime = tr.stats.starttime + (i * fs)
                        val = x[i]
                        duration = count * fs
                        d = {
                            "start_time": starttime.isoformat(),
                            "end_time": endtime.isoformat(),
                            "value": val,
                            "duration": duration,
                            "type": type,
                        }
                        repAmps.append(d)
                        inc = False

            if inc:
                i += 1
    else:
        try:
            import pycheron.repAmpsMetric as repAmpsMetric
        except ImportError:
            subprocess.call(
                [
                    "f2py",
                    "-c",
                    "-m",
                    "repAmpsMetric",
                    os.path.dirname(__file__) + "/repAmps/repAmps.f90",
                ]
            )
            import pycheron.repAmpsMetric as repAmpsMetric

            trunc, start, end, vals, types = repAmpsMetric.repamps(x, minRep, bound)

            for i in range(trunc):
                starttime = tr.stats.starttime + (start[i] * fs)
                endtime = tr.stats.starttime + (end[i] * fs)
                val = vals[i]
                if types[i] == "j":
                    type = "Jitter"
                else:
                    type = "Flat"
                duration = endtime - start
                d = {
                    "start_time": starttime.isoformat(),
                    "end_time": endtime.isoformat(),
                    "value": val,
                    "duration": duration,
                    "type": type,
                }
                repAmps.append(d)

    if generateMasks:
        df = pd.DataFrame(repAmps)
        df.drop("value", 1)
        df.drop("type", 1)
        df.drop("duration", 1)
        m = df.to_dict()
    else:
        m = "No Masks Created"

    out = {
        "snclq": tr.get_id(),
        "count": len(repAmps),
        "rep_amp": repAmps,
        "mask": m,
        "metric_name": "repeatedAmplitudeMetric",
        "start_time": tr.stats.starttime.isoformat(),
        "end_time": tr.stats.endtime.isoformat(),
    }

    return out
