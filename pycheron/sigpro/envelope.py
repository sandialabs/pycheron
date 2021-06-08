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

__all__ = ["envelope"]
import numpy as np
from pycheron.psd.noise.deadChannel import DDT
from scipy.fftpack import hilbert


def envelope(tr):
    """
    Calculates the envelope of a seismic signal for the given input Obpsy Trace

    Python version of envelope in R IRISSeismic package [#]_ , originally adapted from seewave package.
    The envelope is determined by adding the squared amplitudes of the function and it's
    Hilbert-Transform and then taking the square-root. [#]_
    The envelope at the start/end should not be taken too seriously (same process as obspy.signal.filter.envelope), data
    are detrended, demeaned and cosine tapered first. See the DDT function in the psd.noise.deadChannel function for
    more details.

    "The seismic envelope is defined as:

    .. math::

            E(t) = \sqrt{T(t)^2 + H(t)^2}

    where *T(t)* is the seismic trace and *H(t)* is the hilbert transform of *T(t)*" (Callahan, 2020)

    .. note:: Append processing results to trC.stats.processing so can keep track of what's occurred to the trace.

    :param tr: ObsPy Trace object
    :type tr: obspy.core.trace.Trace

    :return: Returns trace whose data have been replaced with the envelope of the seismic signal
    :rtype: obspy.core.trace.Trace

    * Code originally ported from IRISSeismic R Cran Package
      (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
      The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index) and
      augmented and adapted for use within Pycheron.

    **Example**

    .. code-block:: python

        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        from pycheron.sigpro.envelope import envelope

        # Instantiate client object
        client = Client("IRIS")

        # Grab data from 2010-02-27 06:00:00 to 2010-02-27 09:00:00
        t = UTCDateTime("2010-02-27T06:00:00.000")
        st = client.get_waveforms("IU","ANMO","00","BHZ",t,t+180*60)
        tr = st[0]

        #Create trace envelope
        trenv = envelope(tr)

    **Plotting**

    .. code-block:: python

        import matplotlib.pyplot as plt

        # Plot signal and envelope data
        plt.plot(tr.data,color = 'grey', label = 'raw data')
        plt.plot(trenv.data, marker = 'o', markersize = 0.1, linestyle = 'None', label = 'envelope')
        plt.legend()

    .. image:: _static/envelope.png

    .. rubric:: References taken directly from Callahan, 2020

    .. [#] https://cran.r-project.org/web/packages/IRISSeismic/index.html
    .. [#] Kanasewich, E.R., 1981, "Time sequence analysis in geophysics"

    """

    # Copy trace before applying DDT so retain original trace. Detrend, demean, and cosine taper the data. Then take
    # hilbert FFT
    xC = tr.copy()
    xC = DDT(xC)
    hfft = hilbert(xC.data)

    # The envelope is just the modulus of the hfft transformed data
    xC.data = np.sqrt(xC.data ** 2 + hfft ** 2)

    return xC
