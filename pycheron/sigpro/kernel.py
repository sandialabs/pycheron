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

__all__ = ["modified_daniell_kernel"]

import numpy as np
from scipy.fftpack import ifft


def modified_daniell_kernel(m, logger=None):
    """
    Calculates the modified daniell smoothing kernel (discrete symmetric normalized).

    The modified daniell kernel halves the end coefficients. The periodogram is a wildly fluctuating estimate of the
    spectrum with high variance. For a stable estimate, the periodogram must be smoothed. Bloomfield [#]_
    recommends the daniell window as a smoothing filter for generating an estimated spectrum from the periodogram.

    The Daniell filter differs from an evenly weighted moving average (rectangular filter) only in that the first and
    last weights are half as large as the other weights. A plot of the filter weights will thus have the form of a
    trapezoid. The advantage of the Daniell filter over the rectangular filter for smoothing the periodogram is that the
    Daniell filter has less leakage, which refers to the influence of variance at non-Fourier frequencies on the
    spectrum at the Fourier frequencies. The leakage is related to sidelobes in the frequency response of the filter.
    Successive smoothing by Daniell filters with different spans gives an increasingly smooth spectrum, and is
    equivalent to single application of a resultant filter produced by convolution of the individual spans of Daniell
    filters [1]_ [#]_


    .. note::
            Code based on R's ``kernel.R`` code [#]_ . ``kernel.R`` has more kernel options (dirichlet, fejar, daniell),
            but IRIS codes only use the modified daniell smoother so only coded that up with the other needed parts of
            the code.

    :param m: Kernel dimension(s) (e.g., number of weights or span of the filter). When m has a length larger
              than one, it means the convolution of kernels of dimensions ``m[j] for j in 1:length(m)``.
    :type m: int
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger


    :return: returns coefficients of the daniell smoothing kernel
    :rtype: numpy.array

    **Example**

    .. code-block:: python

        #import function
        from pycheron.sigpro.kernel import modified_daniell_kernel

        # Define spans
        spans = 1
        k = modified_daniell_kernel(spans)
        # m = 1
        >>> [0.25, 0.5, 0.25]


        # Define new span where m = 2
        spans = 2
        k = modified_daniell_kernel(spans)
        # m = 2
        >>> [0.125, 0.25, 0.25, 0.25, 0.125]

        # Define array of spans
        spans = np.array([3,5,7,9])
        k = modified_daniell_kernel(spans // 2)
        # m = [1,2,3,4]'
        >>> [0.00016276+7.98571369e-17j 0.00130208+8.57715021e-17j ... 0.00130208+6.14124843e-17j 0.00016276+7.63595486e-17j]


    .. rubric:: References

    .. [#] Bloomfield, P. (1976). Fourier Analysis of Time Series: An Introduction. Wiley, New York.
    .. [#] http://terascan.smast.umassd.edu/nasdata/archive/achaudhu.old/docs/TIMESERIES/notes_6.pdf
    .. [#] https://stat.ethz.ch/R-manual/R-devel/library/stats/html/density.html

    """

    # Check if single number or array
    if isinstance(m, (int, np.int64)):
        # If single number compute coefficients then put kernel coefficients in appropriate order
        k = (np.append(np.repeat(1, m), 0.5)) / (2 * m)
        k = _tskernel(k, m)
        return k
    # If not a single number has to be convolved with kernel, first calculate the modified daniel kernel for first
    # coefficient
    else:
        k = m[0]
        k = modified_daniell_kernel(k)
        # Loop through other coefficients and convolve between input sequence and specific kernel for each of the
        # remaining coefficients
        for i in range(1, len(m)):
            k = _kernapply(
                np.concatenate((np.zeros(m[i]), k, np.zeros(m[i]))),
                modified_daniell_kernel(m[i]),
                m[i],
                logger,
            )
        return k


def _tskernel(k, m):
    """

    Puts kernel coefficients in the appropriate positions

    :param k: initial coefficients
    :param m: kernel dimension(s)

    :return: returns coefficients of the daniel smoothing kernel

    """

    y = []
    # Create descending and ascending array based on m values, e.g., if m = 2, ar =[2,1,0,1,2]
    ar = np.concatenate((range(m, 0, -1), range(0, m + 1)))
    # Loop through ar and assign appropriate weights for each position, e.g., using example above [2,1,0,1,2]
    for i in ar:
        y.append(k[i])
    return y


def _kernapply(x, k, m, logger=None):
    """

    Computes convolution between an input sequence and a specific kernel

    :param x: input vector to be smoothed
    :param k: smoothing object
    :param m: loop value (kernel dimension)

    :return: smoothed version of the input sequence (kernel)

    Note: Took out np.real calculation (originally had this). R code applies real restriction if y is numeric,
    though when circular is used only returns y. As not using non-circular option, simply omitted circular logical input
    argument. Enforcing y to be real is handled in crossSpectrum.py function.

    """

    # Ensure x is longer than kernel
    if len(x) <= (2 * m):
        logger.error("kernapply(): x is shorter than kernel k")
        return
    # If kernel dimension zero, return x
    if m == 0:
        return x
    # Otherwise calculate the weights (w) and the resulting smoothed version of the input sequence (y)
    else:
        n = len(x)
        first = k[m : (2 * m + 1)]
        last = k[0:m]
        w = np.concatenate((first, np.repeat(0, (n - 2 * m - 1)), last))
        y = ifft(np.fft.fft(x) * np.fft.fft(w), n)
        return y
