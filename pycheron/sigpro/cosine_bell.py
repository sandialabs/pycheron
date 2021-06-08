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

__all__ = ["cosine_bell"]

from math import floor
import numpy as np
from pycheron.util.logger import Logger


def cosine_bell(x, taper=0.1, logger=None):
    """
    Apply a cosine bell taper to a time series. Cosine taper is applied to the first and last p[i] observations of time
    series x. Python version of R's spec.taper function that applies a cosine-bell taper to incoming data

    :param x: ObsPy trace.data or any time series
    :type x: numpy.ndarray of floats
    :param taper: The proportion to be tapered at each end of the time series
    :type taper: float
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger

    :return: Tapered time series data in a numpy array
    :rtype: numpy.ndarray

    * Code originally ported from R's spec.taper function
      https://stat.ethz.ch/R-manual/R-devel/library/stats/html/spec.taper.html, which is part of R's
      stats package.  R Core Team (2013). R: A language and environment for statistical computing. R
      Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL
      http://www.R-project.org/. It has been augmented and adapted for use within Pycheron.

    **Example**

    .. code-block:: python

        #Import function
        from pycheron.sigpro.cosine_bell import cosine_bell
        import scipy.signal as sig

        # Initialize IRIS client
        client = Client("IRIS")

        #Define start/end times
        starttime = UTCDateTime("2011-05-01T00:00:00.000")
        endtime = starttime + 3600

        # Grab data from IRIS web client server with specified starttime/endtime and SNCL
        st = client.get_waveforms("CI","PASC","00","BHZ", starttime,endtime)
        print st

        #Grab out trace, define taper
        tr = st[0]
        taper = 0.1

        #Taper data with cosine bell using 0.1 taper
        tD = cosine_bell(tr.data,taper)

        #Normally you would at least detrend then taper. Type = Linear will also demean.
        tD = sig.detrend(tr.data,type='linear')
        tDD = cosine_bell(tr.data,taper)

    **Plotting**

    .. code-block:: python

        import matplotlib.pyplot as plt

        #Plot original data and tapered data
        fig = plt.figure(figsize=(10,10))
        ax1 = fig.add_subplot(311)
        ax1.plot(tr.times(),tr.data)
        ax1.set_title('Original Data')

        ax2 = fig.add_subplot(312)
        ax2.plot(tr.times(),tD)
        ax2.set_title('Detrended-Cosine-tapered data')

        #Plot additional figure
        ax3 = fig.add_subplot(313)
        ax3.plot(tr.times(),tDD)
        ax3.set_title('Cosine-tapered data')

    .. image:: _static/cosine_bell.png

    """
    # TODO: Error handling if length of data == 0

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Ensure taper between 0 and 0.5, otherwise exit
    if taper < 0 or taper > 0.5:
        logger.error("p must be between 0 and 0.5")
        return

    tap = x
    # Get length of time series
    nr = len(x)
    # Get variable length of x * taper
    m = int(floor(nr * taper))
    if m != 0:
        # Calculate taper then apply to time series x
        n = np.arange(1, 2 * m - 1 + 1, 2)
        w = 0.5 - 0.5 * np.cos(np.pi * n / (2 * m))
        a = np.repeat(1, nr - 2 * m)
        ar = np.append(w, a)
        ars = np.append(ar, w[::-1])
        tap = x * ars
    return tap
