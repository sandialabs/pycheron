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

__all__ = ["roll_median"]

import numpy as np


def roll_median(x, nwin=7, increment=1):
    """

    Python version of roll_median function in seismicRoll R package

    Can be utilized to replace outliers detected by the hampel function. Window size is interpeted as the
    full window length. Values within ``n/2`` of the beginning or end of x are set to NA.

    Setting increment to a value greater than one will result in NAs for all skipped over values


    Fast, center-aligned rolling medians algorithm

    :param x: input data
    :type x: numpy.array
    :param nwin: n integer window size
    :type nwin: int
    :param increment: increment integer shift to use when sliding the window to the next location
    :type increment: int

    :return: Returns vector of rolling median values of the same length as x
    :rtype: numpy.array

    **Example**

    .. code-block:: python

        from pycheron.rollseis.roll_median import roll_median
        from pycheron.rollseis.roll_sd import roll_sd
        from psd.noise.findOutliers import findOutliers
        import numpy as np
        # Create sinusoidal signal with outliers
        g =[]
        for i in range(1,101):
            g.append(np.sin(0.1*i))
        # Create indices with outliers to detect
        g=np.asarray(g)
        ginds = [10,19,49,67,68,73,90]
        for i in ginds:
            g[i] = g[i]*10
        g_fixed = np.array(g)
        gm_fixed = np.array(g)

        # Find outliers with rolling_hampel filter
        # Parameters
        nwin = 7
        threshold = 6
        selectivity = None
        increment = 1
        fixedThreshold = True
        outliers = findOutliers(g,nwin,threshold,selectivity,increment,fixedThreshold)

        #Apply rolling median
        for i in outliers:
        g_fixed[i] = roll_median(g,10,1)[i]


    **Plotting**

    .. code-block:: python

        # Using outputs from above example
        import matplotlib.pyplot as plt

        plt.plot(outliers,g[outliers],'ro', label= 'Outliers detected')
        plt.plot(g,color='grey',label='original')
        plt.plot(g_fixed,color = 'red',linestyle = '--',label = 'corrected')
        plt.plot(outliers,g[outliers],'ro')
        plt.ylim(-15, 10)
        plt.title('Outliers detected: %s'%(len(outliers)))
        plt.legend(loc='lower left')

    .. image:: _static/roll_median.png

    """

    # Check if x is not a pd.Series (e.g., Mean from dailyDCOffsetMetric), if not break

    # if not isinstance(x,pd.Series):
    #    print 'the x supplied is not a pd.Series'
    #    exit(1)

    # Check if nwin > length(x) if so break

    if nwin > len(x):
        print("roll_median/medianFilter Error: Nwin cannot be greater than len x ")
        return

    # Avoid infinite loops
    if increment < 1:
        print("roll_median/medianFilter Error: Increment must be >= 1")
        return

    # Get the half-window width
    k = int(nwin // 2)

    # For the valid region, calculate the result
    ind = k

    # Initialize output vector
    out = np.full(len(x), np.nan)

    # For the valid region, calculate the result; e.g., if nwin = 10, k =5, and say x = 100, valid region is between
    # 5 and 95
    while ind < len(x) - k:
        out[ind] = _medianFilter(x, nwin, ind)
        ind += increment
    return out


def _medianFilter(x, nwin, ind):
    """

    Simple rolling median used to replace outliers detected with the hampel filter

    :param x: input data
    :type x: numpy.array
    :param nwin: n integer window size
    :type nwin: int
    :param ind: index looping over
    :type ind: int

    :return: Returns a vector of the same length as incoming data with NaNs where this is not a full window's
             worth of data
    :rtype: numpy.array

    **Example**

    .. code-block:: python


    """

    # Get half-window width
    k = int(nwin // 2)

    # Calculate median value
    tmp = np.full(nwin, np.nan)
    for i in range(nwin):
        idx = int(ind - k + i)
        tmp[i] = x[idx]
    # Sort tmp in correct order
    tmp = sorted(tmp)

    # If nwin even do the following, else, take tmp[nwin/2] to determine outlier replacement
    if nwin % 2 == 0:
        out = (tmp[(nwin // 2) - 1] + tmp[nwin // 2]) / 2
    else:
        out = tmp[nwin // 2]
    return out
