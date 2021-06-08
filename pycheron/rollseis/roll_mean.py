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

__all__ = ["roll_mean"]

import numpy as np
from pycheron.util.logger import Logger


def _meanFilter(x, nwin, ind, align):
    """

    Simple rolling mean with an alignment parameter

    :param x: input data
    :type x: numpy.array
    :param nwin: n integer window size
    :type nwin: int
    :param ind: index looping over
    :type ind: int
    :param align: window alignment, ``"left"``, ``"center"`` or ``"right"``. More detailed information provided in
                  roll_mean function description
    :type align: str

    :return: "Returns a vector of the same length as incoming data with NaNs where this is not a full window's
             worth of data" (Callahan,2020)
    :rtype: numpy.array

    * Code originally ported from seismicRoll R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive R
      Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and
      augmented and adapted for use within Pycheron. This code is equivalent to IRIS's seismicRoll roll_mean in the
      source_files.cpp function

    **Example**

    .. code-block:: python

    """

    # Initialize total to 0
    total = 0

    # Loop through the integer window size
    for i in range(nwin):
        # If index at left edge of window, increment total by ind+integer window size
        if align == "left":
            total += x[ind + i]
        # If index at center of window, halve the window. Increment total by ind - half window + integer window size
        elif align == "center":
            k = nwin // 2
            total += x[ind - k + i]
        # If index at right edge of window, increment total by ind-integer window size
        else:
            total += x[ind - i]
    # Calculate out to be total divided by integer window size
    out = total / nwin
    return out


def roll_mean(x, nwin=7, increment=1, align="center", logger=None):
    """

    Fast rolling mean with alignment parameter.

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
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger

    :return: returns vector of rolling mean values of the same length as x
    :rtype: numpy.array

    .. note:: "Additional performance gains can be achieved by skipping increment values between calculations"
              (Callahan, 2020)

    * Code originally ported from seismicRoll R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll.
      The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and
      augmented and adapted for use within Pycheron. This code is equivalent to IRIS's seismicRoll
      roll_mean_numeric_vector in the source_files.cpp function file. The main roll_mean was skipped because it
      seemed it was really only needed for passing a float to the C code rather than a character code.

    **Example**

    .. code-block:: python

        import numpy as np
        from pycheron.rollseis.roll_mean import roll_mean

        #Create contrived example to show left, right, and center alignment
        x = [1,2,3,4,5]
        x = np.repeat(x,10)

        left = roll_mean(x,5,1,'left')
        right = roll_mean(x,5,1,'right')


    **Plotting**

    .. code-block:: python

        import matplotlib.pyplot as plt

        plt.plot(x, marker = '^', color = 'black', linestyle ='None',label = 'data')

        # Plot left alignment with blue circles using 5 point window
        plt.plot(left, marker = 'o', color = 'blue', markerfacecolor = 'None', label='align=left')

        # Plot right alignment with red circles using 5 point window
        plt.plot(right, marker = 'o', color = 'red', markerfacecolor = 'None', label='align=right')

        # Adjust plot limits to show upper/lower data points more clearly
        plt.xlim([-1,55])
        plt.ylim([0,6])
        #Add legend and title
        plt.legend(loc = 'lower right')
        plt.title('Test of roll_mean with a 5-point window')


    .. image:: _static/roll_mean.png

    """
    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Test if integer window size greater than len x
    if nwin > len(x):
        logger.error("roll_mean(): Nwin cannot be greater than len x")
        return

    # Avoid infinite loops
    if increment < 1:
        logger.error("roll_mean(): Increment must be >= 1")
        return

    # Initialize output vector
    out = np.full(len(x), np.nan)

    # Get the start and endpoints for calculating results in valid region
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
        out[ind] = _meanFilter(x, nwin, ind, align)
        ind += increment
    return out
