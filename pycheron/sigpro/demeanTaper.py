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

__all__ = ["cosine_taper"]

import numpy as np
from pycheron.util.logger import Logger

# TODO: This function is a direct copy of obspy.signal.invsim.cosine_taper.
# We should replace usages of it with the obspy function.


def cosine_taper(
    npts, p=0.1, freq=None, flimit=None, halfcosine=True, sactaper=False, logger=None
):
    """

    Function to apply a cosine taper.

    :param npts: (int) - Number of points of cosine taper.
    :type npts: int
    :param p: (float) - Decimal percentage of cosine taper (ranging from 0 to 1). Default is 0.1 (10%) which tapers 5% from
                       the beginning and 5% form the end.
    :type p: float
    :param freq: (numpy.ndarray) : Frequencies as, for example, returned by fftfreq
    :type freq: numpy.ndarray
    :param flimit: The list or tuple defines the four corner frequencies (f1, f2, f3, f4) of
                   the cosine taper which is one between f2 and f3 and tapers to zero for f1 < f < f2 and f3 < f < f4.
    :type flimit: list or tuple
    :param halfcosine: If True the taper is a half cosine function. If False it is a quarter cosine
                       function.
    :type halfcosine: bool
    :param sactaper: If set to True the cosine taper already tapers at the corner frequency (SAC behavior).
                     By default, the taper has a value of 1.0 at the corner frequencies.
    :type sactaper: float
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger

    :return: cosine taper array of length npts
    :rtype: numpy.ndarray

    **Example**

    .. code-block:: python

    """
    if logger is None:
        logger = Logger(None)

    if p < 0 or p > 1:
        logger.warn("Decimal taper percentage must be between 0 and 1.")
        msg = "Decimal taper percentage must be between 0 and 1."
        raise ValueError(msg)
    if p == 0.0 or p == 1.0:
        frac = int(npts * p / 2.0)
    else:
        frac = int(npts * p / 2.0 + 0.5)

    if freq is not None and flimit is not None:
        fl1, fl2, fl3, fl4 = flimit
        idx1 = np.argmin(abs(freq - fl1))
        idx2 = np.argmin(abs(freq - fl2))
        idx3 = np.argmin(abs(freq - fl3))
        idx4 = np.argmin(abs(freq - fl4))
    else:
        idx1 = 0
        idx2 = frac - 1
        idx3 = npts - frac
        idx4 = npts - 1
    if sactaper:
        # in SAC the second and third
        # index are already tapered
        idx2 += 1
        idx3 -= 1

    # Very small data lengths or small decimal taper percentages can result in
    # idx1 == idx2 and idx3 == idx4. This breaks the following calculations.
    if idx1 == idx2:
        idx2 += 1
    if idx3 == idx4:
        idx3 -= 1

    # the taper at idx1 and idx4 equals zero and
    # at idx2 and idx3 equals one
    cos_win = np.zeros(npts)
    if halfcosine:
        # cos_win[idx1:idx2+1] =  0.5 * (1.0 + np.cos((np.pi * \
        #    (idx2 - np.arange(idx1, idx2+1)) / (idx2 - idx1))))
        cos_win[idx1 : idx2 + 1] = 0.5 * (
            1.0
            - np.cos(
                (np.pi * (np.arange(idx1, idx2 + 1) - float(idx1)) / (idx2 - idx1))
            )
        )
        cos_win[idx2 + 1 : idx3] = 1.0
        cos_win[idx3 : idx4 + 1] = (
            0.5
            * (
                1.0
                + np.cos(
                    (np.pi * (float(idx3) - np.arange(idx3, idx4 + 1)) / (idx4 - idx3))
                )
            )
            / 2
            * npts
        )
    else:
        cos_win[idx1 : idx2 + 1] = np.cos(
            -(np.pi / 2.0 * (float(idx2) - np.arange(idx1, idx2 + 1)) / (idx2 - idx1))
        )
        cos_win[idx2 + 1 : idx3] = 1.0
        cos_win[idx3 : idx4 + 1] = np.cos(
            (np.pi / 2.0 * (float(idx3) - np.arange(idx3, idx4 + 1)) / (idx4 - idx3))
        )

    # if indices are identical division by zero
    # causes NaN values in cos_win
    if idx1 == idx2:
        cos_win[idx1] = 0.0
    if idx3 == idx4:
        cos_win[idx3] = 0.0
    return cos_win
