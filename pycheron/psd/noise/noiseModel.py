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

__all__ = ["noiseModel"]

import numpy as np


def noiseModel(freq):
    """
    The noiseModels function returns the New High Noise Model and New Low Noise Model from the Peterson paper.
    Values are returned for the specific frequencies specified in the freq argument. [#]_ [#]_ [#]_

    New Low/High Noise Model in acceleration (referenced to 1 (m/s^2)^2/Hz):

    .. math::

        NLNM_{acc}/NHNM_{acc} = A + B*\log_{10}(T)

    where, T = minimum period

    If velocity desired:

    .. math::

        NLNM_{vel}/NHNM_{vel} = NLNM_{acc} + 20.0\log_{10}(P/2*\pi)

    Source code to compare:
      * https://github.com/g2e/seizmo/blob/master/noise/nlnm.m
      * https://github.com/g2e/seizmo/blob/master/noise/nhnm.m


    :param freq: array of frequencies
    :type freq: numpy.ndarray

    :return:

        * nlnm - low noise model
        * nhnm - high noise model

    :rtype:

        * numpy.ndarray
        * numpy.ndarray

    **Example**

    .. code-block:: python

        from pycheron.psd.noise.noiseModel import noiseModel
        from pycheron.psd.noise.getNoise import getNoise
        from pycheron.psd.psdList import psdList

        #test data
        data = 'test/test_data/6e_sp06_ehe.407438.tar.mseed'

        # reading in stream
        st = obspy.read(data)

        # calculating psds
        psds = psdList(st)
        # Get instrument corrected psds
        f,n,psd = getNoise(psds)

        freq = f[0]
        period = 1/freq
        # calculating nhnm, nlnm
        nlnm,nhnm = noiseModel(freq)

    **Plotting**

    .. code-block:: python

        import matplotlib.pyplot as plt

        # Normally shown on log xscale need to update
        plt.plot(nlnm[::-1])
        plt.plot(nhnm[::-1])
        plt.title('Noise Models')
        plt.xlabel('Period (s)')
        plt.ylabel('Power (dB)')


    .. image:: _static/noiseModel.png

    .. rubric:: References

    .. [#] Peterson, J, 1993, Observations and Modeling of Seismic Background Noise, U.S.G.S. OFR-93-322
    .. [#] McNamara and Buland, 2003, Ambient Noise Levels in the Continental United States
           (https://pdfs.semanticscholar.org/ed3a/a907fd7a541c8bfb691ddf896df495c406dc.pdf)
    .. [#] McNamara and Boaz, 2005, Seismic Noise Analysis System Using Power Spectral
           Density Probability Density Functions: A Stand-Alone Software Package
           (https://pubs.usgs.gov/of/2005/1438/pdf/OFR-1438.pdf)

    """

    # Obtain period from frequency input
    if isinstance(freq, list):
        freq = np.asarray(freq)
    period = 1 / freq

    # ------------NLMN---------------------------------------------------------
    # Create NLNM table information: minimum period, a, b based on above equation
    minPeriod = np.array(
        [
            0.10,
            0.17,
            0.40,
            0.80,
            1.24,
            2.40,
            4.30,
            5.00,
            6.00,
            10.00,
            12.00,
            15.60,
            21.90,
            31.60,
            45.00,
            70.00,
            101.00,
            154.00,
            328.00,
            600.00,
            10000,
            100000,
        ]
    )
    a = np.array(
        [
            -162.36,
            -166.70,
            -170.00,
            -166.40,
            -168.60,
            -159.98,
            -141.10,
            -71.36,
            -97.26,
            -132.18,
            -205.27,
            -37.65,
            -114.37,
            -160.58,
            -187.50,
            -216.47,
            -185.00,
            -168.34,
            -217.43,
            -258.28,
            -346.88,
            -346.88,
        ]
    )

    b = np.array(
        [
            5.64,
            0.00,
            -8.30,
            28.90,
            52.48,
            29.81,
            0.00,
            -99.77,
            -66.49,
            -31.57,
            36.16,
            -104.33,
            -47.10,
            -16.28,
            0.00,
            15.70,
            0.00,
            -7.61,
            11.90,
            26.60,
            48.75,
            48.75,
        ]
    )

    # Create "breaks" (bins) based on the minimum period, to figure out the appropriate a and b for each period. Return
    # indices of the bins to which each value belongs to
    nlnmBreaks = minPeriod
    nlnmRows = np.digitize(period, nlnmBreaks)

    # Create nlnm using equation NLNMacc = A + B*log10(T) referred to 1 (m/s^2)^2/Hz, where T = period, a and b are
    #     # calculated from above
    # Subtract 1 to agree with R code indexing to ensure grabbing the appropriate numbers
    nlnm = a[nlnmRows - 1] + b[nlnmRows - 1] * np.log10(period)

    # If not acceleration use the following equation: NLNMvel/NHNMvel = NLNMacc + 20.0log10(P/2*pi) db ref 1 (m/sec)^2

    # ------------NHMN---------------------------------------------------------
    # Create NHNM table information: minimum period, a, b based on above equation
    minPeriod = np.array(
        [0.10, 0.22, 0.32, 0.80, 3.80, 4.60, 6.30, 7.90, 15.40, 20.00, 354.80, 100000]
    )

    a = np.array(
        [
            -108.73,
            -150.34,
            -122.31,
            -116.85,
            -108.48,
            -74.66,
            0.66,
            -93.37,
            73.54,
            -151.52,
            -206.66,
            -206.66,
        ]
    )

    b = np.array(
        [
            -17.23,
            -80.50,
            -23.87,
            32.51,
            18.08,
            -32.95,
            -127.18,
            -22.42,
            -162.98,
            10.01,
            31.63,
            31.64,
        ]
    )

    # Create "breaks" (bins) based on the minimum period, to figure out the appropriate a and b for each period. Return
    # indices of the bins to which each value belongs to
    nhnmBreaks = minPeriod
    nhnmRows = np.digitize(period, nhnmBreaks)

    # Create nhnm using equation NHNMacc = A + B*log10(T) referred to 1 (m/s^2)^2/Hz, where T = period, a and b are
    # calculated from above
    # Subtract 1 to agree with R code indexing to ensure grabbing the appropriate numbers
    nhnm = a[nhnmRows - 1] + b[nhnmRows - 1] * np.log10(period)

    return nlnm, nhnm
