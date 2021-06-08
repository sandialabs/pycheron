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

__all__ = ["roll_range"]

import numpy as np
from pycheron.util.logger import Logger


def roll_range(x, nwin=7, increment=1, align="center", logger=None):
    """

    Fast rolling range with an alignment parameter. Python version of roll_range function in main seismicRoll R package.

    :param x: Input data to apply rolling range
    :type x: numpy.array
    :param nwin: nwin Integer window size, interpreted as the full window length.
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

    :return: Returns a vector of rolling range values (difference between max/min values), of the same length as
             incoming data with NaNs where this is not a full window's worth of data.
    :rtype: numpy.array

     .. note:: "Additional performance gains can be achieved by skipping increment values between calculations"
              (Callahan, 2020)

    * Code originally ported from seismicRoll R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive R
      Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and
      augmented and adapted for use within Pycheron. This code is equivalent to IRIS's seismicRoll
      roll_range in the main package.
    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Ensure numpy array has been passed in
    if not isinstance(x, np.ndarray):
        logger.error("roll_range error: x must be a numpy array")

    # Test if integer window size greater than len x
    if nwin > len(x):
        logger.error("roll_range error: Nwin cannot be greater than len x")
        return

    # Avoid infinite loops
    if increment < 1:
        logger.error("roll_range error: Increment must be >= 1")
        return

    # Assign alignment codes for input into _rangeFilter
    if align == "left":
        alignCode = -1
    elif align == "center":
        alignCode = 0
        # Guarantee that nwin is odd
        if nwin % 2 == 0:
            nwin = nwin + 1
    elif align == "right":
        alignCode = 1
    else:
        logger.error("roll_range error: align must be left|center|right")
        return

    # Calculate rolling range
    out = _rangeFilter(x, nwin, increment, alignCode)

    return out


def _rangeFilter(x, nwin, increment, alignCode):
    """
    Simple rolling range with an alignment parameter. Python version of roll_range_numeric_vector function in
    seismicRoll R package.

    :param x: Input data to apply rolling range
    :type x: numpy.array
    :param nwin: nwin Integer window size, interpreted as the full window length.
    :type nwin: int
    :param increment: "Integer shift to use when sliding the window to the next location. Setting increment to a
                             value greater than one will result in NaNs for all skipped over indices." (Callahan, 2020)
    :type increment: int
    :param alignCode: Window alignment code ``"left"``, ``"center"``  or ``"right"`` used to determine start and
                      endpoints for calculating results in valid region
    :type alignCode: str
    :return
    :rtype

    .. note:: "Additional performance gains can be achieved by skipping increment values between calculations"
              (Callahan, 2020)

    * Code originally ported from seismicRoll R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive R
      Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and
      augmented and adapted for use within Pycheron. This code is equivalent to IRIS's seismicRoll
      roll_range_numeric_vector in the source_files.cpp function.
    """

    # Initialize output vector
    out = np.full(len(x), np.nan)

    # Get the start and endpoints for calculating results in valid region
    # "left" aligned window case: index at left edge of window. Start = 0, end = len(x) - (integer window size - 1)
    if alignCode == -1:
        start = 0
        end = len(x) - (nwin - 1)
    # "center" aligned window case: index at center of window. Start = half-width of the integer window size, end =
    # len(x) - (half-width of the integer size window)
    elif alignCode == 0:
        start = nwin // 2
        end = len(x) - nwin // 2
    # "right" aligned window case: index at right edge of window. Start = (integer window size - 1), end = len(x)
    else:
        start = nwin - 1
        end = len(x)

    # Initialize ind at start after set above
    ind = start
    # For the valid region, calculate the rolling range result
    while ind < end:
        out[ind] = _rolling(x, nwin, ind, alignCode)
        ind += increment
    return out


def _rolling(x, nwin, ind, alignCode):
    """
    Calculate range within a rolling window. Python version of roll_range function in source_files.cpp seismicRoll R
    package

    :param x: input data
    :type x: numpy.array
    :param nwin: n integer window size
    :type nwin: int
    :param ind: index looping over
    :type ind: int
    :param alignCode: window alignment, ``"left"``, ``"center"`` or ``"right"``. More detailed information provided in
           roll_range function description
    :type align: str

    :return: Returns a vector of the same length as incoming data with NaNs where this is not a full window's
             worth of data
    :rtype: numpy.array

    * Code originally ported from seismicRoll R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive R
      Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and
      augmented and adapted for use within Pycheron. This code is equivalent to IRIS's seismicRoll
      roll_range in the source_files.cpp function.
    """

    # Set default xMin, xMax
    xMin = x[ind]
    xMax = x[ind]

    # Loop through integer window size
    for i in range(nwin):

        # "left" aligned window case: index at left edge of window. Calculate xMin, xMax
        if alignCode == -1:
            if np.isnan(xMin):
                xMin = x[ind + i]
            if x[ind + i] < xMin:
                xMin = x[ind + i]
            if np.isnan(xMax):
                xMax = x[ind + i]
            if x[ind + i] > xMax:
                xMax = x[ind + i]

        # "center" aligned window case: index at center of window. Calculate xMin, xMax
        elif alignCode == 0:
            k = nwin // 2
            if np.isnan(xMin):
                xMin = x[ind - k + i]
            if x[ind - k + i] < xMin:
                xMin = x[ind - k + i]
            if np.isnan(xMax):
                xMax = x[ind - k + i]
            if x[ind - k + i] > xMax:
                xMax = x[ind - k + i]

        # "right" aligned window case: index at right edge of window. Calculate xMin, xMax
        else:
            if np.isnan(xMin):
                xMin = x[ind - i]
            if x[ind - i] < xMin:
                xMin = x[ind - i]
            if np.isnan(xMax):
                xMax = x[ind - i]
            if x[ind - i] > xMax:
                xMax = x[ind - i]
    # Calculate range value
    out = xMax - xMin

    return out
