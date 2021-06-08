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

__all__ = ["getPDF"]

import numpy as np


def getPDF(noise, lo=-200, hi=-49, binSize=1):
    """
    Computes a probability density function (PDF) as defined by McNamara and Boaz [#]_ from a user defined noise matrix.
    This is Pycheron's version of IRISSeismic's noiseMatrix2PdfMatrix function.

    "The McNamara and Boaz paper describes creating histograms of the discretized power levels at each frequency bin
    associated with a set of PSDs. The value in each cell of the PDF matrix is the fraction of the corrected PSDs that
    have a power level at that frequency bin." (https://cran.r-project.org/web/packages/IRISSeismic/IRISSeismic.pdf)

    .. note::
        Use the default settings to create a PDF matrix that matches the McNamara paper
        (https://cran.r-project.org/web/packages/IRISSeismic/IRISSeismic.pdf)

    :param noise: Noise matrix returned from getNoise function
    :type noise: numpy.ndarray
    :param lo: Lowest frequency bin (power level in dB) for the PDF y-axis (default = -200)
    :type lo: int
    :param hi: Highest frequency bin (power level in dB) for the PDF y-axis. ( **Note:** Python does not include the
     last integer, so if the highest freq you want is -50, you must input -49) (default = -49)
    :type hi: int
    :param binSize: Size in dB of each bin (default = 1)
    :type binSize: int

    :return: Array of probability density values; rows=dB level, columns=frequencies, e.g.,
                        row1 corresponds to lowest frequency bin (e.g., if using default values -200, next row would be
                        -199, etc.)
    :rtype: numpy.ndarray

    * Code originally ported from IRISSeismic R Cran Package
       (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
       The Comprehensive R Archive Network.
       Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index) and
       augmented and adapted for use within Pycheron

    **Example**

    .. code-block:: python

        from pycheron.psd.psdList import psdList
        from pycheron.psd.psdStatistics import psdStatistics
        from pycheron.psd.McNamaraPSD import McNamaraPSD
        from pycheron.psd.McNamaraBins import McNamaraBins
        from pycheron.psd.getPDF import getPDF
        from pycheron.psd.noise.getNoise import getNoise
        from matplotlib.ticker import ScalarFormatter
        from obspy.imaging.cm import pqlx

        #test data
        data = 'test/test_data/7a_cabn_bhe.884965.tar.mseed'

        #reading in stream
        st = obspy.read(data)

        tr=st[0]

        # creating psds
        psd=psdList(st)
        # calculating psd statistics
        stats = psdStatistics(psd)

        # getting noise, frequency
        f,n,psd = getNoise(psd)

        # function parameters
        lo = -200 #Default
        hi = -51 #Default
        binSize = 1 #Default
        freq = n[0]
        pdf = getPDF(freq, lo, hi, binSize)



    **Plotting**

    .. code-block:: python

        # using outputs from above example
        period=1/stats[0]['noise_matrix_frequency']
        pdf_rev=np.fliplr(pdf)
        fig1, ax1 = plt.subplots()
        plt.pcolor(pdf_rev,cmap=pqlx,alpha=1.0,vmin=0,vmax=30)
        plt.xscale('log')
        ax1.get_xaxis().set_major_formatter(ScalarFormatter())
        ax1.grid(True,which='both')
        ax1.set_ylabel('Amplitude [$m^2/s^4/Hz$] [dB]')
        ax1.set_xlabel('Period (s)')
        plt.colorbar()
        plt.show()

    .. image:: _static/getPDF.png

    .. rubric:: References

    .. [#] McNamara and Boaz, 2005, Seismic Noise Analysis System Using Power Spectral
           Density Probability Density Functions: A Stand-Alone Software Package
           (https://pubs.usgs.gov/of/2005/1438/pdf/OFR-1438.pdf)
    """

    # Create array of bins from hi to lo; add one to hi bin because np.arange is not endpoint inclusive.
    bins = np.arange(lo, hi + 1, binSize)
    nbins = list(range(len(bins) - 1))

    # Binning noise array -- digitize returns the indices of the bins to which each value in the input array belongs.
    # So discrete will contain the index value of the bin it occupies at that frequency for that particular PSD segment.
    discrete = np.digitize(noise, bins)

    # Initialize numpy array container for pdfMatrix based on the length of nbins and the columns in discrete (# of
    # frequencies). Total length should be number of bins (nbins = bins -1) x number of frequencies
    pdfMatrix = np.empty((len(nbins), discrete.shape[1]))

    # Looping through columns of discrete
    # Make discrete into a list so we can count the number of hits per bin
    for i in range(discrete.shape[1]):
        cols = discrete[:, i].tolist()

        # looping through bins, counting the number of hits per bin. Put the count value into the appropriate bin and
        # frequency of the pdfMatrix
        for j in nbins:
            c = cols.count(nbins[j])
            pdfMatrix[j][i] = c

    return pdfMatrix
