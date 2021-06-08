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

__all__ = ["roll_sd"]

import numpy as np
from math import sqrt


def _sdFilter(x, nwin, ind, align):
    """

    Simple rolling standard deviation with an alignment parameter. Python version of roll_sd function in
    source_files.cpp seismicRoll R package

    :param x: (vector) - input data

    :param nwin: (int) - n integer window size, interpreted as the full window length

    :param ind: (int) - index looping over

    :param align: (str) "window alignment, "left", "center" (default), or "right". Note for align = "center" the window
                       size is increased by one if necessary to guarantee an odd window size. The align parameter
                       determines the alignment of the current index within the window. Thus:

                        align = "left"  [*------] will cause the returned vector to have n-1 NaN values at the right end
                        align = "center" [---*---] will cause the returned vector to have (n-1)/2 NaN values at either
                                end
                        align = "right" [------*] will cause the returned vector to have n-1 NaN values at the left end"
                        (Callahan, 2020)

    :return: (vector) - out - A vector of rolling standard deviation values of the same length as x

    * Code originally ported from seismicRoll R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive R
      Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and
      augmented and adapted for use within Pycheron. This code is equivalent to IRIS's seismicRoll
      roll_sd in the source_files.cpp function.

    """

    # Initialize total to 0
    total = 0
    # Loop through integer window size
    for i in range(nwin):
        # Index at left edge of window, increment total by ind+integer window size
        if align == "left":
            total += x[ind + i]
        # Index at center of window, halve the window. Increment total by ind - half window + integer window size
        elif align == "center":
            k = nwin // 2
            total += x[ind - k + i]
        # Index at right edge of window, increment total by ind - integer window size
        else:
            total += x[ind - i]

    # Calculate mean
    mean = total / nwin
    # Initialize out
    out = 0
    # Loop through integer window size
    for i in range(nwin):
        # Index at left edge of window
        if align == "left":
            out += (x[ind + i] - mean) * (x[ind + i] - mean)
        # Index at center of window
        elif align == "center":
            k = nwin // 2
            out += (x[ind - k + i] - mean) * (x[ind - k + i] - mean)
        # Index at right edge of window
        else:
            out += (x[ind - i] - mean) * (x[ind - i] - mean)
    # Get std
    out = sqrt(out / (nwin - 1))

    return out


def roll_sd(x, nwin, increment, align):
    """

    Fast rolling standard deviation with an alignment parameter. Python version of roll_sd_numeric_vector function in
    source_files.cpp seismicRoll R package

    :param x: Input data
    :type x: numpy.array
    :param nwin: n Integer window size, interpreted as the full window length.
    :type nwin: int
    :param increment: "Integer shift to use when sliding the window to the next location. Setting increment to a
                             value greater than one will result in NaNs for all skipped over indices." (Callahan, 2020)
    :type increment: int
    :param align: "Window alignment, ``"left"``, ``"center"`` (default), or ``"right"``. **Note:** for
                  ``align = "center"`` the window size is increased by one if necessary to guarantee an odd window
                  size. The align parameter determines the alignment of the current index within the window. Thus:

                  * ``align = "left"``  [\*------] will cause the returned vector to have n-1 NaN values at the right
                    end
                  * ``align = "center"`` [---\*---] will cause the returned vector to have ``(n-1)/2`` NaN values at
                    either end
                  * ``align = "right"`` [------\*] will cause the returned vector to have n-1 NaN values at the left
                    end"
                  (Callahan, 2020)

    :type align: str

    :return: returns vector of rolling standard deviation values of the same length as x
    :rtype: numpy.array

    * Code originally ported from seismicRoll R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive R
      Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and
      augmented and adapted for use within Pycheron. This code is equivalent to IRIS's seismicRoll
      roll_sd_numeric_vector in the source_files.cpp function.

    **Example**


    """
    # Initialize output vector with nans; make it the length of x
    out = np.full(len(x), np.nan)

    # Index at left edge of window. Start = 0, end = len(x) - (integer window size - 1)
    if align == "left":
        start = 0
        end = len(x) - (nwin - 1)
    # Index at center of window. Start = half-width of the integer window size,
    # end = len(x) - (half-width of the integer size window)
    elif align == "center":
        start = nwin // 2
        end = len(x) - nwin // 2
    # Index at right edge of window. Start = integer window size -1. End = len(x)
    else:
        start = nwin - 1
        end = len(x)

    # # Initialize ind at start after set above
    ind = start

    # For the valid region, calculate the rolling std result
    while ind < end:
        out[ind] = _sdFilter(x, nwin, ind, align)
        ind += increment
    return out
