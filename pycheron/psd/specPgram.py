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

import numpy as np
import pylab
from numba import njit

__all__ = ["specPgram"]


def specPgram(x, fs, taper=0.1, fast=True, demean=True, detrend=True):
    """
    Python version of R's spec.pgram that calculates a Periodigram using a FFT, though with more limited input options.

    :param x: trace data
    :type x: numpy.ndarray
    :param fs: sampling frequency
    :type fs: float
    :param taper: taper
    :type taper: float
    :param fast: fast fourier transfer using highly composite number.
    :type fast: bool
    :param demean: demean the trace
    :type demean: bool
    :param detrend: detrend the data
    :type detrend: bool

    :return: an array of frequencies and an array of spectral values
    :rtype: numpy.ndarray

    * Code originally ported from R's spec.pgram function
      https://stat.ethz.ch/R-manual/R-devel/library/stats/html/spec.pgram.html, which is part of R's
      stats package.  R Core Team (2013). R: A language and environment for statistical computing. R
      Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL
      http://www.R-project.org/.It has been augmented and adapted for use within Pycheron. For example, the Python
      version does not include the span and kernel input options.

    **Example**

    .. code-block:: python

        from pycheron.psd.specPgram import specPgram
        import obspy

        #test data
        data = 'test/test_data/7a_cabn_bhe.884965.tar.mseed'

        #reading in stream
        st = obspy.read(data)

        #getting trace
        tr = st[0]

        # getting trace data values
        x = tr.data

        #sampling rate
        fs = tr.stats.sampling_rate

        #taper
        taper=.1

        #FFT
        fast = True

        #demean
        demean=True

        #detrend
        detrend=True

        freq, pxx = specPgram(x, fs, taper, fast, demean, detrend)

        print 'frequency:', freq #Frequency
        >>> frequency: [1.15620303e-05 2.31240606e-05 3.46860909e-05 ... 1.99999769e+01 1.99999884e+01 2.00000000e+01]
        print 'pxx:', pxx #spectral values
        >>> pxx: [9.82917567e+09 1.02135447e+10 6.66068275e+09 ... 8.93294708e-04 1.65414423e-04 5.15396272e-04]

    """
    # Assign sampling freq
    xfreq = fs
    # Convert x to an array
    x = np.asarray(x)
    # Get length of x
    n = len(x)

    # detrend
    if detrend:
        x = _specDetrend(x)

    # demean
    if demean:
        x = pylab.detrend_mean(x, axis=0)

    # taper
    x = _specTaper(x, taper)

    # Correction for tapering
    u2 = 1 - (5.0 / 8.0) * taper * 2

    # If fast = True, then get the next highly compositable number for computing fft
    if fast:
        n = _nextN(n)

    # Return floor of next highly compositable number. Divide freq by that to get step increment for calculating new
    # frequency array
    nspec = np.floor(n / 2)
    step = float(xfreq) / n
    # Calculate frequency using npsec and step
    freq = np.arange(start=step, stop=(nspec * step) + step, step=step)
    # Get fft of x
    xfft = np.fft.fft(x)
    # Get periodogram
    pgram = xfft * np.conj(xfft) / (n * xfreq)

    pgram = pgram[1 : int(nspec + 1)]

    # Grab out the real values and correct for tapering
    spec = pgram.real
    spec = spec / u2

    return freq, spec


@njit(nogil=True)
def _specDetrend(x):
    """
    Python version of R's detrend function found in spec.pgram. Detrends data

    :param x: numpy array
    :return: detrended array

     * Code originally ported from R's spec.pgram function
       https://stat.ethz.ch/R-manual/R-devel/library/stats/html/spec.pgram.html, which is part of R's
       stats package.  R Core Team (2013). R: A language and environment for statistical computing. R
       Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL
       http://www.R-project.org/. It has been augmented and adapted for use within Pycheron.

    """
    # Get length of x
    n = len(x)
    t = np.asarray(list(range(1, n + 1))) - (n + 1) / 2.0

    sumt2 = n * (n ** 2 - 1) / 12.0
    # Detrend
    d = (x - np.mean(x)) - (np.sum(x * t) * (t / sumt2))
    return d


def _specTaper(x, p=0.1):
    """

    Python version of R's spec.taper function that applies a cosine-bell taper to incoming data

    :param x: numpy array
    :param p: percent to taper (0.0 - 0.5) each end of the incoming data

    :return: taper

    * Code originally ported from R's spec.taper function
      https://stat.ethz.ch/R-manual/R-devel/library/stats/html/spec.taper.html, which is part of R's
      stats package.  R Core Team (2013). R: A language and environment for statistical computing. R
      Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL
      http://www.R-project.org/. It has been augmented and adapted for use within Pycheron.

    """

    # Get len of x
    n = len(x)
    # Calculate floor of len x times percent to taper
    m = np.floor(n * p)
    # Get cosine taper window
    w = 0.5 * (1 - np.cos(np.pi * np.asarray(list(range(1, int(2 * m), 2))) / (2 * m)))
    c = np.append(w, np.ones(int(n - 2 * m)))
    c = np.append(c, np.flip(w, 0))
    # Taper x with the cosine bell taper
    x = c * x
    return x


def _nextN(n):
    """

    Python version of R's nextn function that finds the next highly compositable number for computing fft.

    :param n: length of array

    :return: the smallest integer, greater than or equal to n,
             which can be obtained as a product of powers of the values in the factors list. Default factors are
             hard-coded to [2, 3, 5]

    * Code originally ported from R's nextn function
      https://stat.ethz.ch/R-manual/R-devel/library/stats/html/nextn.html, which is part of R's
      stats package.  R Core Team (2013). R: A language and environment for statistical computing. R
      Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL
      http://www.R-project.org/. It has been augmented and adapted for use within Pycheron.

    """
    # Set up empty values list
    vals = []
    # Set up default set of postivie integer factors
    factors = [2, 3, 5]
    # Loop through factors and while v is less than the length of the array find the smallest integer, greater than or
    # equal to n that can be obtained as a product of powers of the values in the factors list
    for f in factors:
        count = 0
        v = 0
        while v <= n:
            v = count ** f
            count += 1
        vals.append(v)
    return min(vals)
