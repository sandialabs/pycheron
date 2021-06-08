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

__all__ = ["roll_stalta"]
import os
import numpy as np
import subprocess
from pycheron.util.logger import Logger

from pathlib2 import Path

# Example on how to recomplie seismicRoll lib
# ----------------------
# from numpy import f2py
# with open("/Users/jbobeck/pycheron/rollseis/seismicRoll.f", "r") as myfile:
#     source = myfile.read()
# f2py.compile(source, "seismicRoll")


import platform

if platform.system() != "Windows":
    try:
        import pycheron.seismicRoll as seismicRoll
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
                "seismicRoll",
                "--quiet",
                os.path.dirname(__file__) + "/seismicRoll.f90",
            ]
        )
        os.chdir(pwd)
        try:
            import pycheron.seismicRoll as seismicRoll
        except ImportError:
            import seismicRoll


def roll_stalta(x, n_sta, n_lta, increment=1, fortran=False, logger=None):
    """

    Simple rolling STA/LTA ratio calculation utilized for automatic detection of seismic signal arrival times. Python
    version of roll_stalta_numeric_vector function in source_files.cpp seismicRoll R package

    "``roll_stalta`` doesn't do any preprocessing of incoming data and merely calculates the ratio of the average value
    in the STA window to the average value in the LTA window. Windows are aligned so that the index is at the left edge
    of the STA window and at the right edge of the LTA window, e.g., [#]_

    .. math::

        STA(x_{i}) = (1/ns)* \sum_{j = i}^{i+ns} x_{i}

        LTA(x_{i}) = (1/nl) * \sum_{j=i-nl}^{i} x_{i}

        r_{i} = STA_{i}/LTA_{i}

    .. code-block:: console

        [---------- LTA --------*]........
        .......................[*- STA --]

    For proper use of this algorithm seismic data should be preprocessed in following manner:

    demean, detrend and taper the raw signal
    square the processed signal to get power" (Callahan, 2020)

    :param x: input data vector
    :type x: numpy.array
    :param n_sta: integer STA window size
    :type n_sta: int
    :param n_lta: integer LTA window size
    :type n_lta: int
    :param increment: "increment shift to use when sliding the window to the next location. For increments greater
                      than one, the rolling means will not align properly, hence the need for a dedicated
                      ``roll_stalta`` function. Setting increment to a value greater than 1 will result in NaNs for
                      all skipped over indices." (Callahan, 2020)
    :type increment: int
    :param fortran: Whether to use Fortran or not. **Note:** Linux/iOS only
    :type fortran: bool
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger

    :return: Returns vector of values of the same length as x with each point containing the STA/LTA ratio at that point
    :rtype: numpy.array

    .. note:: "Values within n_lta - 1 of the beginning and n_sta - 1 of the end are set to NaNs" (Callahan, 2020)

    * Code originally ported from seismicRoll R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive R
      Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and
      augmented and adapted for use within Pycheron. This code is equivalent to IRIS's seismicRoll
      roll_stalta_numeric_vector in the source_files.cpp function.

    **Example**

    .. code-block:: python

        import numpy as np
        from pycheron.rollseis.roll_stalta import roll_stalta
        import obspy
        from pycheron.psd.noise.deadChannel import DDT

        #Contrived example:
        x = [1,5,3,2,1]
        x = np.repeat(x,20)

        #calculate rolling_stalta with n_sta = 3, n_lta = 6, increment = 1 for above vector x
        p = roll_stalta(x,3,6)

    **Plotting**

    .. code-block:: python

        import matplotlib.pyplot as plt
        # Plot original data as black triangles
        plt.plot(x,marker='^',color = 'black', linestyle = 'None', label = 'data')

        #Plot rolling_stalta in red circles
        plt.plot(p, marker = 'o', color = 'red', markerfacecolor = 'None', label = 'STA/LTA')

        #Adjust plot limits to show upper/lower data points more clearly
        plt.xlim([-1,101])
        plt.ylim([0,6])
        #Add title and legend
        plt.title('Test of roll_stalta on contrived example')
        plt.legend(loc = 'upper right')

    .. image:: _static/roll_stalta.png

    .. rubric:: References

    .. [#]  http://en.wikipedia.org/wiki/First_break_picking
    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # If using fortran
    if fortran:

        # Ensure integer window sizes not greater than len x for n_sta and n_lta
        if n_sta > len(x):
            logger.error("roll_stalta(): n_sta cannot be greater than len x")
            return
        elif n_lta > len(x):
            logger.error("roll_stalta(): n_lta cannot be greater than len x")
            return

        # Avoid infinite loops
        if increment < 1:
            logger.error("roll_stalta(): Increment must be >= 1")
            return

        # Initialize output vector and call the fortran version of the code base to calculate the rolling STA/LTA
        out = seismicRoll.roll_stalta(x, n_sta, n_lta, increment)
        out[out == 0] = np.nan

    # Otherwise still ensure the same thing but don't need to call the fortran version of the code base
    else:
        # Ensure integer window sizes not greater than len x
        if n_sta > len(x):
            logger.error("roll_stalta(): n_sta cannot be greater than len x")
            return
        elif n_lta > len(x):
            logger.error("roll_stalta(): n_lta cannot be greater than len x")
            return

        # Avoid infinite loops
        if increment < 1:
            logger.error("roll_stalta(): Increment must be >= 1")
            return

        # Initialize output vector to the length of x and fill with nans
        out = np.full(len(x), np.nan)

        # Set ind to n_lta
        ind = n_lta

        # For valid region, calculate the rolling STA/LTA result
        while ind < (len(x) - n_sta):
            out[ind] = stalta_python(x, n_sta, n_lta, ind)
            ind += increment

    return out


def stalta_python(x, n_sta, n_lta, ind):
    """

    This is a simple ratio of two rolling means that is used to detect seismic signal arrivals. Python
    version of roll_stalta function in source_files.cpp seismicRoll R package

    IN:

    :param x: (vector) - input data

    :param n_sta: (int) - integer STA window size (sec)

    :param n_lta: (int) - integer LTA window size (sec)

    :param ind: (int) - index looping over

    OUT:

    :return: (vector) - out - "Returns a vector of the same length as the incoming data with NaNs in the LTA window
                              length at the left end and in the STA window length at the right end" (Callahan, 2020)

    * Code originally ported from seismicRoll R Cran Package
      (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive R
      Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index) and
      augmented and adapted for use within Pycheron. This code is equivalent to IRIS's seismicRoll
      roll_stalta in the source_files.cpp function.

    """

    # Initialize total to 0
    total = 0

    # Loop through integer STA window size and update the total to index + integer STA window size.
    # Calculate the STA, aligned so that the window is right of ind (including ind). STA is total / integer STA window
    # size
    for i in range(n_sta):
        total += x[ind + i]
    sta = total / n_sta

    # Initialize total again for LTA
    total = 0

    # Loop through integer STA window size and update the total to index - integer LTA window size.
    # Calculate the LTA, aligned so that the window is left of ind (including ind). LTA is total/ integer LTA window
    # size
    for i in range(n_lta):
        total += x[ind - i]
    lta = total / n_lta

    # Calculate STA/LTA ratio
    out = sta / lta
    return out
