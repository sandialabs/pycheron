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

__all__ = ["McNamaraPSD"]

from math import floor, log
import numpy as np
from pycheron.psd.McNamaraBins import McNamaraBins
from pycheron.psd.specPgram import specPgram


def McNamaraPSD(tr, loFreq=0.005, hiFreq=10, alignFreq=0.1, binned=True):
    """

    Computes the Power Spectral Density (PSD) for the provided trace, using McNamara's algorithm in [#]_.

    :param tr: ObsPy Trace object
    :type tr: obspy.core.trace.Trace
    :param loFreq: Low end of frequency binning range
    :type loFreq: int
    :param hiFreq: High end of frequency binning range
    :type hiFreq: int
    :param alignFreq: Alignment frequency for determining frequency bins
    :type alignFreq: int
    :param binned: Boolean to determine whether to bin data or not
    :type binned: bool

    :return: list containing the following outputs:

             * freq (`numpy.ndarray`) - frequency
             * spec (`numpy.ndarray`)- spectrum
             * snclq (`list`) - ["station.network.channel.location.quality"]
             * starttime (`str`) - Start time
             * endtime (`str`) - Endtime

    :rtype: list

    * Code originally ported from IRISSeismic R Cran Package
      (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
      The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index) and
      augmented and adapted for use within Pycheron

    **Methodology**

    "The PSD algorithm is designed to be used on one to three hour segments of seismic data and will return output
    containing the (potentially binned) spectrum for that segment. The McNamara algorithm is similar to MATLAB's
    pwelch [#]_ function and has the following steps:

    1. Calculate the average spectrum

       a. Truncate incoming segment of trace data to the nearest power of 2 samples.

       b. Divide each truncated segment into 13 chunks with 75% overlap
          First chunk begins at 0/16 and ends at 4/16. The 13'th chunk begins at 12/16 and ends at 16/16.
          The chunks overlap in the following manner

          .. code-block:: console

                1---5---9---3---
                 2---6---0---
                  3---7---1---
                   4---8---2---

       c. Demean, detrend, and taper the chunk

       d. Calculate the one-sided spectrum for the chunk

       e. Average all 13 spectra to get an averaged spectrum

    2. Create smoothed version of spectrum with binning

       When ``binned = True`` , McNamara style binning is turned on and a smoothed spectrum is returned that contains
       many fewer points than the full spectrum. When these arguments are not specified, binning is automatically turned
       off and the full spectrum is returned.

       Frequencies for binning are generated at 1/8 octave intervals aligned to alignFreq. The power (dB) associated
       with each frequency bin is calculated by averaging over the entire octave centered at that frequency.

       **Note:** The spectra returned by ``McNamaraPSD()`` have not had instrument correction applied.

    3. Convert binned spectra to decibels" (https://cran.r-project.org/web/packages/IRISSeismic/IRISSeismic.pdf)

    **Example**

    .. code-block:: python

        #import function
        import matplotlib.pyplot as plt
        from pycheron.psd.McNamaraPSD import McNamaraPSD

        #test data
        data = 'test/test_data/7a_cabn_bhe.884965.tar.mseed'

        #reading in stream
        st = obspy.read(data)

        #getting trace
        tr = st[0]

        #defining parameters
        loFreq =.005 #Default
        hiFreq =10 #Default
        alignFreq = .1 #Defualt
        binned = True #Defualt

        psd = McNamaraPSD(tr, loFreq=.005, hiFreq=10, alignFreq=.1, binned=True)

        #Flags to obtain this data are (note these are uncorrected!!!!):
        #   psd[0] = frequency
        #   psd[1] = psd
        #   psd[2] = snclq
        #   psd[3] = Starttime of data segment
        #   psd[4] = Starttime of data segment

        # Output not shown for brevity
        print 'Frequency:', psd[0] #Frequency
        print 'PSD:', psd[1] #PSDs
        print 'SNCLQ:', psd[2] #SNCLQ
        print 'Starttime:', psd[3] #Starttime of data segment
        print 'Endtime:', psd[4] #Endtime of data segment

    .. rubric:: References

    .. [#] McNamara, D. E., & Boaz, R. I. (2006). Seismic noise analysis system using power spectral density probability
           density functions: A stand-alone software package. US Geological Survey.
    .. [#] https://www.mathworks.com/help/signal/ref/pwelch.html

    """
    # Get data, start, end times as well as general snclq information
    data = tr.data
    starttime = tr.stats.starttime
    endtime = tr.stats.endtime
    # if file is not mseed, it will not have data quality flag, but first try to get it. If fail, just get the sncl
    try:
        snclq = [
            "{0}.{1}.{2}.{3}.{4}".format(
                str(tr.stats.network),
                str(tr.stats.station),
                str(tr.stats.location),
                str(tr.stats.channel),
                str(tr.stats.mseed.dataquality),
            )
        ]
    except AttributeError:
        snclq = [
            "{0}.{1}.{2}.{3}".format(
                str(tr.stats.network),
                str(tr.stats.station),
                str(tr.stats.location),
                str(tr.stats.channel),
            )
        ]
    # Truncate trace data to the nearest power of 2 samples
    pow2 = floor(log(len(data), 2))
    truncateLen = 2 ** pow2

    # Initialize container, get trace sampling rate
    specSum = np.zeros(1)
    fs = tr.stats.sampling_rate

    # Divide each truncated trace into 13 segments with 75% overlap
    for i in range(13):
        first = int(i * (truncateLen / 16))
        last = int((i + 4) * truncateLen / 16)

        d = data[first:last]

        # checking if len() of d is prime, if so -1
        # This is done because of the FFT done in the periodogram.
        if len(d) % 2 == 1:
            d = d[:-1]

        # Calculate Periodogram
        f, pxx = specPgram(d, fs=fs)
        # Scale periodogram with 10% taper scale factor
        specSum = specSum + 2 * pxx

    # average the Spectral sum
    pxx = specSum / 13

    # Bin if binned=True
    if binned:
        freq, spec = McNamaraBins(f, pxx, loFreq, hiFreq, alignFreq)
    # Otherwise return raw spectrum
    else:
        freq = f
        spec = pxx

    # Convert to dB
    spec = 10 * np.log10(spec)

    # Return output list
    psd = [freq, spec, snclq, starttime.isoformat(), endtime.isoformat()]

    return psd
