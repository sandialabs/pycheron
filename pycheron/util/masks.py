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

__all__ = ["consecutive", "samples2time"]


def consecutive(data, stepsize=1):
    """
    Cluster consecutive elements where indices differ by stepsize. Utility to group masks/QC issue start/stop times

    :param data: array of indices to cluster consecutive elements for trace
    :type data: numpy.ndarray
    :param stepsize: stepsize to use for clustering consecutive elements where indices differ.
    :type stepsize: int

    :return: Returns list of arrays of consecutive indices
    :rtype: list
    """

    return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)


def samples2time(mask_index, fs, start):
    """
    This function takes in index locations and converts the locations into time

    :param mask_index: array containing the index locations of where there is a QC issue
    :type mask_index: numpy.array
    :param fs: Sampling frequency
    :type fs: float
    :param start: Start time of tr
    :type start: obspy.core.utcdatetime.UTCDateTime

    :return: List of dictionaries with the following keys and types:

            * start_time (`str`) - Start time of mask
            * end_time (`str`) - End time of mask

    :rtype: list
    """

    # Set index to mask_index provided
    index = mask_index
    # Calculate offsets based on index * 1/sampling rate
    offsets = mask_index * 1.0 / fs

    # Initialize start and out list
    st = start
    out = []

    # loop through the offsets and grab out start and end times for mask indices
    for i in range(len(offsets) - 1):
        starttime = (st + offsets[i]) - (1.0 / fs)
        if index[i] + 1 != index[i + 1]:
            endtime = st + offsets[i]
            out.append({"start_time": starttime.isoformat(), "end_time": endtime.isoformat()})
        elif i != 0 and index[i] - 1 != index[i - 1]:
            inc = 0
            while i + inc < len(index) and index[i] + inc == index[i + inc]:
                inc = inc + 1
                endtime = st + offsets[i + (inc - 1)]
                out.append(
                    {
                        "start_time": starttime.isoformat(),
                        "end_time": endtime.isoformat(),
                    }
                )
        elif i == 0:
            inc = 0
            while i + inc < len(index) and index[i] + inc == index[i + inc]:
                inc = inc + 1
                endtime = st + offsets[i + (inc - 1)]
                out.append(
                    {
                        "start_time": starttime.isoformat(),
                        "end_time": endtime.isoformat(),
                    }
                )

    return out
