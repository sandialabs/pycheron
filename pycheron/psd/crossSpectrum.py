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

__all__ = ["crossSpectrum"]

import numpy as np
import scipy.signal as sig
from math import floor
from obspy.core.trace import Trace
from pycheron.sigpro.cosine_bell import cosine_bell
from pycheron.sigpro.kernel import modified_daniell_kernel, _kernapply


def crossSpectrum(
    x, y, fs, spans=None, taper=0.1, pad=0, demean=False, detrend=True, logger=None
):
    """
    CrossSpectrum is a function specifically for calculating cross-spectral values between two time series.

    :param x: Time series trace objecy
    :type x: obspy.core.trace.Trace
    :param y: Time series trace object
    :type y: obspy.core.trace.Trace
    :param fs: Sampling frequency of the x and y time series
    :type fs: int
    :param spans: Array of odd integers that provide widths of modified Daniell smoothers to be used to smooth
                  the periodogram (i.e. np.array([3,5]) <- default span provided in transferFunctionMetric)
    :type spans: numpy.ndarray
    :param taper: Specifies the proportion of data to taper. A split cosine bell taper is applied to this
                  proportion of the data at the beginning and end of the series. Percentage in decimal format
                  (i.e. .1 = 10%)
    :type taper: float
    :param pad: Proportion of data to pad. Zeros are added to the end of the series to increase its length by the
                proportion pad
    :type pad: int
    :param demean: If True, subtract the mean of the series (DEFAULT = False)
    :type demean: bool
    :param detrend: If True, remove a linear trend from the series. (Default = True). Linear will also demean.
    :type detrend: bool
    :param logger: Logger object
    :type logger: pycheron.util.logger.Logger

    :return: Returns the following parameters and types:

        * freq (`numpy.ndarray`) - spectral frequencies
        * spec1 (`numpy.ndarray`) - 'two-sided' spectral amplitudes for x (real-valued only)
        * spec2 (`numpy.ndarray`) - 'two-sided' spectral amplitudes for y (real-valued only)
        * coh (`numpy.ndarray`)- magnitude squared coherence between x and y
        * phase (`numpy.ndarray`) - cross-spectral phase between x and y
        * Pxx (`numpy.ndarray`) - periodogram for x
        * Pyy (`numpy.ndarray`) - periodogram for y
        * Pxy (`numpy.ndarray`) - cross-periodogram for x & y

    .. note::

        This function marries functionality from R's spec.pgram [#]_ [#]_  [#]_ with MATLAB's pwelch [#]_. Other
        features mimic the code in Octave's pwelch() function [#]_.

    **Example**

    .. code-block:: python

        # Import function
        from pycheron.psd.crossSpectrum import crossSpectrum
        from math import pi
        import matplotlib.pyplot as plt

        # Initialize IRIS client
        client = Client("IRIS")

        #Define start/end times
        starttime = UTCDateTime("2011-05-01T00:00:00.000")
        endtime = starttime + 3600

        # Grab data from IRIS web client server with specified starttime/endtime and SNCL
        st1 = client.get_waveforms("CI","PASC","00","BHZ", starttime,endtime)
        st2 = client.get_waveforms("CI","PASC","10","BHZ", starttime, endtime)

        # Grab out traces
        tr1 = st1[0]
        tr2 = st2[0]

        sampling_rate = int(round(tr1.stats.sampling_rate))

        # Calculate the cross spectrum
        freq, spec1, spec2, coh, phase, pxx, pyy, pxy = crossSpectrum(tr1, tr2, sampling_rate, spans=np.array([3,5,7,9]), taper=.1, pad=0, demean=False, detrend=True)

        # Calculate the transfer function
        transferFunction = pxy/pxx
        transferAmp = abs(transferFunction)
        transferPhase = phase * 180 /pi

    **Plotting**

    .. code-block:: python

        # Plot transfer function amplitude and phase from above example
        fig = plt.figure(figsize=(10,10))
        ax1 = fig.add_subplot(211)
        ax1.semilogx(1/freq,transferAmp)
        ax1.set_title('Transfer Function Amplitude')
        ax1.set_xlabel('Period (sec)')

        ax2 = fig.add_subplot(212)
        ax2.semilogx(1/freq,transferPhase)
        ax2.set_title('Transfer Function Phase')
        ax2.set_xlabel('Period (sec)')
        ax2.set_ylabel('Degrees')

    .. image:: _static/crossSpectrum.png

    .. rubric:: References

    .. [#] https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/spectrum
    .. [#] https://cran.r-project.org/web/packages/psd/vignettes/normalization.pdf
    .. [#] http://www.ltrr.arizona.edu/~dmeko/notes_6.pdf
    .. [#] https://www.mathworks.com/help/signal/ref/pwelch.html
    .. [#] http://octave-signal.sourcearchive.com/documentation/1.0.7/pwelch_8m-source.html
    .. [#] http://terascan.smast.umassd.edu/nasdata/archive/achaudhu.old/docs/TIMESERIES/notes_6.pdf


    """

    # If both spans and kernel are None, notify the user they must specify at least one valid span or kernel

    if spans is None:
        logger.error(
            "crossSpectrum(): Must specify spans in order for smoothed periodogram estimate to be calculated"
        )
        return

    m = spans // 2

    # TODO: Handle Error if spans size is zero.

    if spans.size != 0:
        kernel = modified_daniell_kernel(m)

    # check lengths of data to make sure they are equal
    n = len(x)
    nspec = int(floor(n / 2))
    xfreq = fs
    freq = np.linspace(xfreq / n, (xfreq / n) * nspec, nspec)

    # need to get data out of traces x and y so that taper doesn't error
    # if detrend and demean are false. cosine_bell takes a numpy array. detrend
    # and demean functions convert obspy trace to numpy array of data.
    # If detrend and demean are not run, need to have a numpy array not a trace.
    if all(isinstance(inst, Trace) for inst in [x, y]):
        x = x.copy().data
        y = y.copy().data

    # TODO: The DDT function in psd/noise/deadChannel.py seems to implement the
    # functionality to demean, detrend, and taper traces. Could we use it instead
    # of the following three if statements?
    if detrend:
        x = sig.detrend(x, type="linear")
        y = sig.detrend(y, type="linear")

    if demean:
        x = sig.detrend(x, type="constant")
        y = sig.detrend(y, type="constant")

    # TODO: handle negative or greater than 0.5 values in taper parameter.
    # negative values will return data currently, but should raise exception.
    # greater than 0.5 values will throw exception in fast fourier transform,
    # but probably should throw exception here so it is more clear.
    # # taper
    if taper > 0:
        x = cosine_bell(x, taper, logger)
        y = cosine_bell(y, taper, logger)

    # to correct for tapering: Bloomfield (1976, p. 194). Total taper is taper * 2
    u2 = 1 - (5 / 8) * taper * 2

    # If need to pad time series, create array of zeros based on n * pad, then pad to end of the time series.
    # Recalculate n, nspec, and frequencies based on new array length
    if pad > 0:
        xa = np.zeros(n * pad)
        x = np.append(x, xa)
        y = np.append(y, xa)
        n = len(x)
        nspec = int(floor(n / 2))
        freq = np.linspace(xfreq / n, (xfreq / n) * nspec, nspec)

    # Calculate x periodogram
    fft_x = np.fft.fft(x)
    pxx = fft_x * np.conj(fft_x) / (n * xfreq)
    pxx[0] = 0.5 * (pxx[1] + pxx[n - 1])

    # y periodogram
    fft_y = np.fft.fft(y)
    pyy = fft_y * np.conj(fft_y) / (n * xfreq)
    pyy[0] = 0.5 * (pyy[1] + pyy[n - 1])

    # cross spectrum
    pxy = fft_x * np.conj(fft_y) / (n * xfreq)
    pxy[0] = 0.5 * (pxy[1] + pxy[n - 1])  # TODO: is this correct?

    # Apply kernel to smooth periodogram
    if kernel.size != 0:
        pxx = _kernapply(pxx, kernel, sum(m), logger)
        pyy = _kernapply(pyy, kernel, sum(m), logger)
        pxy = _kernapply(pxy, kernel, sum(m), logger)

    # force spectra to be real (as took out np.real in kernel calculation, see kernel.py for comments)
    pxx = pxx[1 : int(nspec + 1)]
    spec1 = np.real(pxx[0 : int(nspec)])
    pyy = pyy[1 : int(nspec + 1)]
    spec2 = np.real(pyy[0 : int(nspec)])
    pxy = pxy[1 : int(nspec + 1)]

    # Calculate coherence and phase
    coh = abs(pxy) ** 2 / (spec1 * spec2)
    phase = np.angle(pxy)

    # Apply correction for tapering
    spec1 /= u2
    spec2 /= u2

    return freq, spec1, spec2, coh, phase, pxx, pyy, pxy
