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

# Code base originally written by Steven Magana-Zook from Lawrence Livermore National Laboratories
# Augmented and adapted for use within Pycheron by Pycheron team

__all__ = ["find_max_discontinuity"]

import numpy as np


def find_max_discontinuity(trace, factor=10, max_range=1e5):
    """

    Function to find the max discontinuity within a given trace object

    :param tr: obspy trace object
    :type tr: obspy.core trace object
    :param factor: scale factor used to scale the calculated differences standard deviation (default = 10) by to create
                   a threshold for determining extreme values. This value is used as the discontinuity threshold.
    :type factor: int
    :param max_range: Maximum range for the trace's data
    :type max_range: float

    return: For each extreme value, the standard deviation from the rest of the trace is calculated. This returns the
            max standard deviation of the extreme values found in the trace.

    rtype: float

    * Code base originally written by Steven Magana-Zook from Lawrence Livermore National Laboratories
      Augmented and adapted for use within Pycheron by Pycheron team

    """

    # Calculate range of data
    mrange = np.max(trace.data) - np.min(trace.data)

    # If mrange is greater than max range, scale the trace data by max_range/mrange
    if mrange > max_range:
        scale = max_range / mrange
        trace.data *= scale

    # Subtract each i+1 element from each ith element then square those differences
    diffs = np.zeros(len(trace.data))
    differences = diffs[0 : len(diffs) - 1] = (
        trace.data[1:] - trace.data[0 : len(trace.data) - 1]
    )
    differences = differences ** 2

    # Compute statistics on the differences to generate a discontinuity threshold
    diff_avg = np.mean(differences)
    diff_std = np.std(differences)
    threshold = factor * diff_std

    # Call any value an "extreme value" if it's value minus the mean is above the threshold
    extreme_values = differences[differences - diff_avg > threshold]

    # Find out how many standard deviations away each extreme value is from the rest of the trace
    std_away = np.array(
        [0 if diff_std == 0 else num / diff_std for num in extreme_values]
    )

    # Return the max standard deviation
    return np.max(std_away)
