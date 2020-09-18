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

    Simple rolling standard deviation with an alignment parameter

    :param x: (vector) - input data

    :param nwin: (int) - n integer window size

    :param ind: (int) - index looping over

    :param align: (str) window alignment, "left", "center" (default), or "right". Note for align = "center" the window
                       size is increased by one if necessary to guarantee an odd window size. The align parameter
                       determines the alignment of the current index within the window. Thus:

                        align = "left"  [*------] will cause the returned vector to have n-1 NaN values at the right end
                        align = "center" [---*---] will cause the returned vector to have (n-1)/2 NaN values at either
                                end
                        align = "right" [------*] will cause the returned vector to have n-1 NaN values at the left end

    :return: (vector) - out - A vector of rolling standard deviation values of the same length as x

    """

    # Calculate the mean (including the ind)
    total = 0
    for i in range(nwin):
        # Index at left edge of window
        if align == "left":
            total += x[ind + i]
        # Index at center of window
        elif align == "center":
            k = nwin // 2
            total += x[ind - k + i]
        # Index at right edge of window
        else:
            total += x[ind - i]

    mean = total / nwin
    out = 0
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

    out = sqrt(out / (nwin - 1))

    return out


def roll_sd(x, nwin, increment, align):
    """

    Fast rolling standard deviation with an alignment parameter

    :param x: Input data
    :type x: numpy.array
    :param nwin: n Integer window size. The window size nwin is interpreted as the full window length.
    :type nwin: int
    :param increment: Integer shift to use when sliding the window to the next location. Setting increment to a
                             value greater than one will result in NaNs for all skipped over indices.
    :type increment: int
    :param align: Window alignment, ``"left"``, ``"center"`` (default), or ``"right"``. **Note:** for
                  ``align = "center"`` the window size is increased by one if necessary to guarantee an odd window
                  size. The align parameter determines the alignment of the current index within the window. Thus:

                  * ``align = "left"``  [\*------] will cause the returned vector to have n-1 NaN values at the right
                    end
                  * ``align = "center"`` [---\*---] will cause the returned vector to have ``(n-1)/2`` NaN values at
                    either end
                  * ``align = "right"`` [------\*] will cause the returned vector to have n-1 NaN values at the left end

    :type align: str

    :return: returns vector of rolling standard deviation values of the same length as x
    :rtype: numpy.array

    **Example**

    .. code-block:: python

    """
    # Initialize output vector
    out = np.full(len(x), np.nan)

    # Index at left edge of window
    if align == "left":
        start = 0
        end = len(x) - (nwin - 1)
    # Index at center of window
    elif align == "center":
        start = nwin // 2
        end = len(x) - nwin // 2
    # Index at right edge of window
    else:
        start = nwin - 1
        end = len(x)

    ind = start
    # For the valid region, calculate the result
    while ind < end:
        out[ind] = _sdFilter(x, nwin, ind, align)
        ind += increment
    return out
