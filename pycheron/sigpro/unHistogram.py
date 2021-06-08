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

__all__ = ["unHistogram"]

import numpy as np


def unHistogram(vec, startVal=1, inc=1):
    """

    Calculates associated bin values with the proper count of each value when provided with a histogram array or
    ordered set of binned counts of incrementing values (ascending)

    :param vec: A histogram array or ordered set of binned counts
    :type vec: numpy.array
    :param startVal: Initial value of the first bin element (DEFAULT = 1)
    :type startVal: int
    :param inc: Increment rate of each subsequent bin
    :type inc: int

    :return: Bin values with appropriate counts of each
    :rtype: numpy.array

    * Code originally ported from IRISSeismic R Cran Package
      (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
      The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index) and
      augmented and adapted for use within Pycheron.

    **Example**

    .. code-block:: python


    """

    # Initialize array container
    k = np.array([])

    # Use ndenumerate, multidimensional index iterator to loop through vec and append elements to array to get bin
    # values with the appropriate counts of each
    for ind, i in np.ndenumerate(vec):
        if i > 0:
            k = np.append(k, (np.repeat(startVal + (ind[0] * inc), i)))
    if len(k) % 10 != 0:
        k = np.append(k, k[-1])
    try:
        k.reshape(np.int(np.round(len(k) / 10)), 10)
    except ValueError:
        while len(k) != 470:
            k = np.append(k, k[-1])

        k.reshape(np.int(np.round(len(k) / 10)), 10)

    return k
