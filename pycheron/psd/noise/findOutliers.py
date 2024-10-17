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

__all__ = ["hampel", "hampelFilter", "findOutliers"]

import numpy as np
import os
import subprocess
from pathlib2 import Path
import platform
from pycheron.util.logger import Logger

# import hampelF fortran code to speed up processing
if platform.system() != "Windows":
    try:
        import pycheron.hampel as hampelF
    except ImportError:
        fpath = os.path.dirname(__file__)
        pwd = os.getcwd()
        if pwd != str(Path(fpath).parent.parent):
            os.chdir(str(Path(fpath).parent.parent))
        print("------------------------------------------------------")
        print("Building Fortran Library")
        print("------------------------------------------------------")
        subprocess.call(["f2py", "-c", "-m", "hampel", "--quiet", fpath + "/outliers/outliers.f"])
        os.chdir(pwd)
        try:
            import pycheron.hampel as hampelF
        except ImportError:
            import hampel as hampelF


def hampel(data, nwin, ind):
    """
    Python Version of the roll_hampel function in seismicRoll R package

    :param data: input data
    :type data: numpy.ndarray
    :param nwin: n integer window size
    :type nwin: int
    :param ind: index looping through
    :type ind: int

    :return: 1D array of values of the same length as data
    :rtype: numpy.ndarray

    * Code originally ported from IRISSeismic R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive
      R Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and augmented and
      adapted for use within Pycheron

    **Example**

    .. code-block:: python

    """

    # Define the Mean Absolute Deviation Constant
    L = 1.4826

    # Obtain half-window width
    k = nwin // 2

    # Calculating x0 = median(data[nwin])
    # Set tmp variable based on index - half window width to index - half window width + integer window size
    tmp = data[ind - k : ind - k + nwin]
    tmp = np.sort(tmp)

    # If window modulo 2 is 0
    if nwin % 2 == 0:
        x0 = (tmp[(nwin // 2) - 1] + tmp[nwin // 2]) // 2
        # Calculate the absolute difference between each point and the median
        absMinusMedian = abs(data[ind - k : ind - k + nwin] - x0)

        # Calculate the median of the two values after sorting
        absMinusMedian = np.sort(absMinusMedian)
        medianAbsMinusMedian = (absMinusMedian[(nwin // 2) - 1] + absMinusMedian[nwin // 2]) // 2
    # Otherwise
    else:
        x0 = tmp[nwin // 2]
        # Calculate the absolute difference between each point and the median
        absMinusMedian = abs(data[ind - k : ind - k + nwin] - x0)

        # Calculate the median of the two values after sorting
        absMinusMedian = np.sort(absMinusMedian)
        medianAbsMinusMedian = absMinusMedian[nwin // 2]

    # Calculate S0 which is L * the median of teh two values
    S0 = L * medianAbsMinusMedian

    # Create output
    out = abs(data[ind] - x0) / S0
    return out


def hampelFilter(data, nwin, increment):
    """
    Python version of roll_hampel_numeric_vector in seismicRoll R package.
    "The Hampel filter is a robust outlier detector using Median Absolute Deviation (MAD)"
    (https://cran.r-project.org/web/packages/seismicRoll/seismicRoll.pdf); vector of values returned
    that can be tested against different threshold values. "Higher values in return are associated with a higher
    likelihood that the associated point is an outlier when compared with it's neighbors"
    (https://cran.r-project.org/web/packages/seismicRoll/seismicRoll.pdf).

    :param data: input data
    :type data: numpy.ndarray
    :param nwin: n integer window size
    :type nwin: int
    :param increment: Increment integer shift to use when sliding the nwin to the next location
    :type increment: int

    :return: Returns a vector/array of same length of incoming data with NANs in the half-window at either end
    :rtype: numpy.ndarray

    * Code originally ported from IRISSeismic R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive
      R Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and augmented and
      adapted for use within Pycheron

    .. warning::

        "The default value of **increment** = 1 should not be changed. Outliers are defined as individual points that
        stand apart from their neighbors. Applying the Hampel filter to every other point by using
        increment > 1 will invariably miss some of the outliers."
        https://cran.r-project.org/web/packages/seismicRoll/seismicRoll.pdf



    """

    # initialize output vector with NaNs
    d = np.full(len(data), np.nan)

    # Obtain half-window width
    k = int(nwin // 2)

    # Calculate result
    for ind in range(k, len(data) - k, increment):
        d[ind] = hampel(data, nwin, ind)
    return d


def findOutliers(
    data,
    nwin=41,
    threshold=10,
    selectivity=None,
    increment=1,
    fixedThreshold=True,
    fortran=False,
    logger=None,
):
    """
    "Outlier detection with rolling hampel filter; A wrapper for the roll_hampel() function that counts outliers using
    either a user specified threshold value or a threshold value based on the statistics of the incoming data"
    (https://cran.r-project.org/web/packages/seismicRoll/seismicRoll.pdf)


    :param data: (vector/array) - input data
    :type data: numpy.ndarray
    :param nwin: (int) - n integer window size
    :type nwin: int
    :param threshold: (int) minimum value for outlier detection.
    :type threshold: int
    :param selectivity: (float) - value between [0-1] used in determining outliers, or None if fixedThreshold = True.
    :type selectivity: float
    :param increment: (int) - increment integer shift to use when sliding the nwin to the next location.
    :type increment: int
    :param fixedThreshold: (boolean) - Boolean specifying whether outlier detection uses selectivity
    :type fixedThreshold: bool
    :param fortran: (boolean) - Boolean specifying whether to use Fortran libs or not. If libs will not compile or on a
                                Windows Machine, set to False
    :type fortran: bool

    :return: returns a list of indices associated with outliers in the incoming data
    :rtype: list

    * Code originally ported from SeismicRoll R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive
      R Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and augmented and
      adapted for use within Pycheron

    .. note::

          "The threshold and selectivity parameters work like squelch and volume on a CB radio: threshold sets a noise
          threshold below which you don't want anything while selectivity increases the number of points defined as
          outliers. Of course nwin, the window Size, is important as well.


          For B* and V* channels, the default value threshold = 6.0 seems to work well, while for L* channels a value of
          threshold = 12 seems more reasonable. More testing is needed."
          (https://cran.r-project.org/web/packages/seismicRoll/seismicRoll.pdf)

    **Parameter Notes**

    "* threshold (`int`)

      * Threshold level is similar to a sigma value for normally distributed data. Hampel filter values
        above 6.0 indicate a data value that is extremely unlikely to be part of a normal distribution
        (~1/500 million) and therefore likely to be an outlier. By choosing a relatively large value for
        threshold min one can make it less likely that we will generate false positives.
        False positives can include high frequency environmental noise.

    * selectivity (`float`)

      * The selectivity is a value between 0 and 1 and is used to generate an appropriate threshold for
        outlier detection based on the statistics of the incoming data. A lower value for selectivity will
        result in more outliers while a value closer to 1.0 will result in fewer.

    * increment `(int`)

      * The default value of increment=1 should not be changed. Applying the Hampel filter to every other
        point by using increment > 1 will invariably miss some of the outliers."
        (https://cran.r-project.org/web/packages/seismicRoll/seismicRoll.pdf)

    **Example**

    .. code-block:: python

        #import function
        from pycheron.psd.noise.findOutliers import findOutliers

        #Create sinusoidal signal with outliers
        g =[]
        for i in range(1,101):
            g.append(np.sin(0.1*i))
        #Create indices with outliers to detect
        g=np.asarray(g)
        ginds = [10,19,46,67,68,73,90]
        for i in ginds:
            g[i] = g[i]*10


        #Parameters
        nwin = 7
        threshold = 6
        selectivity = None
        increment = 1
        fixedThreshold = True

        #finding outliers
        outliers = findOutliers(g,nwin,threshold,selectivity,increment,fixedThreshold)

    **Plotting**

    .. code-block:: python

        import matplotlib.pyplot as plt

        plt.plot(g,color='grey')
        plt.plot(outliers,g[outliers],'ro')
        plt.ylim(-15, 10)
        plt.title('Outliers detected: %s'%(len(outliers)))

    .. image:: _static/outliersPlot.png

    """
    # Set up logger
    if logger is None:
        logger = Logger(None)

    # If using fortran, call the fortran hample filter lib, set nans if needed
    if fortran:

        h = hampelF.hampelfilter(data, nwin, increment)
        h[h == -99999] = np.nan
        h[h == np.inf] = np.nan

    # Otherwise use the python version to calculate
    else:
        h = hampelFilter(data, nwin, increment)

        # Set to nan if H blows up to infinity. Will happen if 50% of values in window are same
        h[h == np.inf] = np.nan

    # If all nans, exit finding outliers
    if np.isnan(h).all():
        logger.error("roll hampel returns a vector with all NaNs; 50% of values in all windows identical")
        return

    # Find maxH in the array
    maxH = np.nanmax(h)
    # Check if values cross threshold
    if maxH < threshold:
        # If no values cross threshold, return an empty vector
        h = []
        return h

    # Otherwise
    else:
        # Some values cross threshold, set up a new threshold based on maxH and selectivity
        # If fixedThreshold set to True, pick outliers based on threshold, otherwise pick outliers based on maxH and
        # selectivity
        if fixedThreshold is True:
            outinds = [inds for inds, value in enumerate(h) if value > threshold]
            # mask = h > threshold
            # outliers = h[mask]
        else:
            # Added in returning indices in case important for use in plotting, other , etc
            outinds = [inds for inds, value in enumerate(h) if value > maxH * selectivity]
            # mask = h > (maxH * selectivity)
            # outliers = h[mask]

        return outinds
