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

# -*- coding: utf-8 -*-


__all__ = ["psdStatistics"]

from pycheron.psd.noise.getNoise import getNoise
from pycheron.psd.noise.noiseModel import noiseModel
from pycheron.psd.noise.noiseModelInfra import noiseModelInfra
from pycheron.psd.getPDF import getPDF
import numpy as np
from pycheron.util.logger import Logger


def psdStatistics(psds, evalresp=None, logger=None, database=None, special_handling=None):
    """
    Calculates basic statistics (e.g., mean, median, max) as well as envelopes (95,5; 90, 10) of a set of PSDs.

    :param psds: list of psds, use ``psdList()`` function output as input
    :type psds: list
    :param evalresp: numpy.ndarray of freq, amp, phase information matching output of
        client.evalresp, optional as will grab the response information from IRIS if
        not supplied by user. If supplied by user, user should adhere to client.evalresp
        output format
    :type evalresp: numpy.ndarray
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger
    :param database: Database object
    :type database: pycheron.db.sqllite_db.Database
    :param special_handling: Switches on customized handling for data other than seismic. 
        Can be either None or 'infrasound'. None option should be used for seismic data, and 'infrasound' option
        for infrasound data. Other options could be added in the future. Infrasound option could be used for 
        hydrophone data as well (for response removal), but also would compare PSD PDFs to IDC infra models.
    :type special_handling: string

    :return: List of dictionaries (one for each trace in stream) containing the
             following keys and types:

             * snclq (`list`) - list of SNCLQ information, i.e. ``["s.n.c.q.l"]``
             * metric_name (`string`) - name of metric
             * start_time (`string`) - psd start time
             * end_time (`string`) = psd end time
             * noise_matrix_Noise (`numpy.ndarray`): Noise – array of corrected power levels
             * noise_matrix_frequency (`numpy.ndarray`): Frequency – array of associated frequencies for noise array
             * pdf_matrix (`numpy.ndarray`) – array of probability density values; rows=dB level, columns=frequencies
             * pdf_bins (`numpy.array`) – vector of power values (dB) associated with pdfMatrix rows
             * max (`numpy.array`) – maximum power level at each frequency
             * min (`numpy.array`) – minimum power level at each frequency
             * mean (`numpy.array`) – mean power level at each frequency
             * median (`numpy.array`) – median power level at each frequency
             * mode (`numpy.array`) – mode of power level at each frequency (obtained from pdfMatrix)
             * nlnm (`numpy.array`) – low noise model power level at each frequency. 
               # If special_handling = 'infrasound' uses IDC low noise from Brown et. al. 2012 as described in
               # noiseModelInfra. Otherwise uses Peterson models as described in noiseModel. 
             * nhnm (`numpy.array`) – high noise model power level at each frequency
               # If special_handling = 'infrasound' uses IDC high noise from Brown et. al. 2012 as described in
               # noiseModelInfra. Otherwise uses Peterson models as described in noiseModel. 
             * percent_above_nhnm (`numpy.array`) – percent of PSDs above the high noise model at each frequency
               # If special_handling = 'infrasound' uses IDC high noise from Brown et. al. 2012 as described in
               # noiseModelInfra. Otherwise uses Peterson models as described in noiseModel. 
             * percent_below_nlnm (vector) – percent of PSDS below the low noise model at each frequency
               # If special_handling = 'infrasound' uses IDC low noise from Brown et. al. 2012 as described in
               # noiseModelInfra. Otherwise uses Peterson models as described in noiseModel. 
             * percent_95 (`numpy.array`) - 95 percentile power level at each frequency
             * percent_90 (`numpy.array`) - 90 percentile power level at each frequency
             * percent_10 (`numpy.array`) - 10 percentile power level at each frequency
             * percent_5 (`numpy.array`) - 5 percentile power level at each frequency

    :rtype: list of dict

    * Code originally ported from IRISSeismic R Cran Package
      (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
      The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index) and
      augmented and adapted for use within Pycheron. For example, envelopes added by Pycheron team.

    **Example**

    .. code-block:: python

        # import function
        import obspy
        import matplotlib.pyplot as plt
        from pycheron.psd.psdList import psdList
        from pycheron.psd.psdStatistics import psdStatistics

        # test data
        data1 = 'test/test_data/7a_cabn_bhe.884965.tar.mseed'

        # reading in stream
        st = obspy.read(data1)

        # creating list of psds, this is a REQUIRMENT as the psdStatistic takes in a psdList
        psds = psdList(st)

        # Calculating stats
        stats = psdStatistics(psds)

        # psdStatistics returns a list of dictionaries. Each dictionary corresponds to a trace in the st, if there is
        # only one dict in the list, then there is only one trace in the stream.
        # Data can be accessed using the following dict keys

        # stats[0]['noise_matrix_frequency']
        # stats[0]['noise_matrix_noise']
        # stats[0]['PDFMatrix']
        # stats[0]['PDFBins']
        # stats[0]['Max']
        # stats[0]['Min']
        # stats[0]['Mean']
        # stats[0]['Median']
        # stats[0]['Mode']
        # stats[0]['NHNM']
        # stats[0]['NLNM']
        # stats[0]['Percent_above_NHNM']
        # stats[0]['Percent_below_NLNM']
        # stats[0]['95percent']
        # stats[0]['90percent']
        # stats[0]['10percent']
        # stats[0]['5percent']

        # output not show for brevity
        print 'NoiseMatrix:', stats[0]['noise_matrix_frequency'][0] #Frequency matrix of first segment of trace
        print 'NoiseMatrix:', stats[0]['noise_matrix_noise'][0] #Noise matrix of first segment trace
        print 'PDFMatrix:', stats[0]['PDFMatrix'][0] #PDF matrix of first segment of trace
        print 'PDFBins:', stats[0]['PDFBins'][0] #PDFBins (power) associated with PDFMatrix of first segment of trace
        print 'Max:', stats[0]['Max']  #Maximum power level at each frequency
        print 'Min:', stats[0]['Min']  #Minimum power level at each frequency
        print 'Mean:', stats[0]['Mean']  #Mean power level at each frequency
        print 'Median:', stats[0]['Median']  #Median power level at each frequency
        print 'Mode:', stats[0]['Mode']  #Mode power level at each frequency
        print 'Percent_above_NHNM:', stats[0]['Percent_above_NHNM']  #Percent of PSDs above NHNM/IDC high noise model
        print 'Percent_below_NLNM:', stats[0]['Percent_below_NLNM']  #Percent of PSDs below NLNM/IDC low noise model
        print '95percent:', stats[0]['95percent']  #95 percentile
        print '90percent:', stats[0]['90percent']  #90 percentile
        print '10percent:', stats[0]['10percent']  #10 percentile
        print '5percent:', stats[0]['5percent']  #5 percentile

    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Hack -- change later (maybe?) If receive empty psds exit out with k
    k = []
    # Check if psd input is empty or zero, if so return empty list.
    try:
        len(psds[0])
        if len(psds[0][0]) == 0:
            logger.error("psdStatistics(): psd input is empty or zero")
            return k
    except TypeError:
        logger.error("psdStatistics(): psd input is empty or zero")
        return k
    except IndexError:
        logger.error("psdStatistics(): psd input is empty or zero")
        return k

    # Get list of frequencies for associated PSDs as well as noise matrices for each corrected PSD
    # (for each trace in stream)
    if psds and special_handling is None:
        f, n, psds = getNoise(psds, evalresp, logger, units='acc')
    elif psds and special_handling is 'infrasound':
        f, n, psds = getNoise(psds, evalresp, logger, units='def')
    else:
        logger.error("psdStatistics(): psd input is empty or zero")
        return k

    # Doing this again because getNoise can also return an empty psdList
    if not f:
        logger.error("psdStatistics(): psd input is empty or zero")
        return k

    # Initialize data container
    data = []

    # ------------Calculate PSD Max, Min, Mean, Median, Envelopes ------------------------------------------------

    # loop through psds and grab out appropriate frequencies and noise
    for h in range(len(psds)):
        # Get freq, noise
        freq = f[h]
        noise = n[h]

        # Define rows and columns based on dimensions of noise array
        nrow = noise.shape[0]
        ncol = noise.shape[1]

        # initialize stat and envelope containers with NaNs for storage based on ncol value from noise matrix
        # (ncol = # of freqs)
        colMax = np.empty(ncol) * np.nan
        colMin = np.empty(ncol) * np.nan
        colMean = np.empty(ncol) * np.nan
        colMedian = np.empty(ncol) * np.nan
        col95 = np.empty(ncol) * np.nan
        col05 = np.empty(ncol) * np.nan
        col90 = np.empty(ncol) * np.nan
        col10 = np.empty(ncol) * np.nan

        # Loop through rows of columns and fill in one frequency at a time for PSD max, min, mean, median, and envelopes
        for i in range(ncol):
            colMax[i] = np.nanmax(noise[:, i])
            colMin[i] = np.nanmin(noise[:, i])
            colMean[i] = np.nanmean(noise[:, i])
            colMedian[i] = np.nanmedian(noise[:, i])
            col05[i] = np.nanpercentile(noise[:, i], 5)
            col95[i] = np.nanpercentile(noise[:, i], 95)
            col90[i] = np.nanpercentile(noise[:, i], 90)
            col10[i] = np.nanpercentile(noise[:, i], 10)

        # ------------Calculate pctAboveNHNM and pctBelowNLNM-----------------------
        # Obtain high and low noise model based on specific frequencies requested in freq argument. If special handling is
        # infrasound, use IDC noise models, otherwise use Peterson noise models
        if special_handling is 'infrasound':
            nlnm, nhnm = noiseModelInfra(freq)
        else:
            nlnm, nhnm = noiseModel(freq)

        # "Calculate percent metrics based on frequencies that have a noise model value and noise value -- this is in R
        # code though I'm not sure how often this will actually be different than ncol but added it anyway
        # High and low noise models have the same upper and lower frequency limits, so only check one model. Also only
        # need to check one noise PSD as should be the same length" (Callahan, 2020)
        try:
            notNan = np.where(~np.isnan(noise[1, :]) & ~np.isnan(nlnm))
        except IndexError:
            notNan = np.where(~np.isnan(noise[0]) & ~np.isnan(nlnm))

        # Use notNan to construct the aboveCount/belowCount numpy array containers if not equal to ncol
        if len(notNan[0]) == ncol:
            aboveCount = np.zeros(ncol)
            belowCount = np.zeros(ncol)
        else:
            aboveCount = np.empty(ncol) * np.nan
            aboveCount[notNan] = 0
            belowCount = np.empty(ncol) * np.nan
            belowCount[notNan] = 0

        # Loop through notNan (will be equivalent to rows of columns for majority of cases) and determine which values
        # are above or below nhnm and nlnm. The expressions results in a boolean array with the same shape as noise,
        # sums over those boolean elements treating True values as 1 and False values as 0. Does this for each array of
        # noise to get a full aboveCount and belowCount.
        for j in range(len(notNan[0])):
            aboveCount[j] = np.sum((noise[:, j] > nhnm[j]))
            belowCount[j] = np.sum((noise[:, j] < nlnm[j]))

        # Calculate percent of total values above/below to get the final pctAbove nhnm and pctBelow nlnm values
        pctAbove = 100 * (aboveCount / nrow)
        pctBelow = 100 * (belowCount / nrow)

        # ------------Calculate Mode------------------------------------------------

        # For mode, we need to convert to the discretized pdfMatrix that contains the count of PSDs with a given value
        # within each bin. The mode then is just the bin with the highest count, as the mode is the value that occurs
        # most frequently.
        # Good Default dbBins: lo=-200, hi=-50, binsize=1, which is like McNamara psds for seismic, dbBins for infra might be better set at lo=-100, hi=40, binsize=1
        # Calculate hi and lo values based on noise matrix, use binsize = 1 (hard-coded)
        try:
            lo = np.floor(np.nanmin(noise))
            hi = np.ceil(np.nanmax(noise))

            # just incase lo and hi == nan or inf
            if np.isnan(lo) or np.isinf(lo) and special_handling is None:
                lo = -200
            elif np.isnan(lo) or np.isinf(lo) and special_handling is 'infrasound':
                lo = -100
            if np.isnan(hi) or np.isinf(hi) and special_handling is None:
                hi = -50
            elif np.isnan(hi) or np.isinf(hi) and special_handling is 'infrasound':
                hi = 40

        except ValueError:
            lo = -200
            hi = -50

        # Add one because np.arange is not endpoint inclusive
        pdfBins = np.arange(lo, hi + 1, 1)

        # Take noise matrix and compute a probability density function (PDF) based on McNamara and Boaz method
        pdfMatrix = getPDF(noise, lo - 0.5, hi + 0.5, 1)

        # Initialize containers with zeros for storage based on ncol value from noise matrix (# of frequencies)
        colMode = np.zeros(ncol)

        # Loop through columns to calculate mode
        for k in range(ncol):
            # looping through rows of column and calculate mode for each bin, e.g., bin with highest count
            for m in range(len(pdfMatrix[:, k])):
                if pdfMatrix[:, k][m] == np.nanmax(pdfMatrix[:, k]):
                    colMode[k] = pdfBins[m]

        # ------------Output--------------------------------------------------------
        # Stuff output into dictionary and append to final data list object
        d = {
            "metric_name": "psdStatistics",
            "snclq": psds[h][0][2][0],
            "start_time": psds[h][0][3],
            "end_time": psds[h][-1][4],
            "noise_matrix_frequency": freq,
            "noise_matrix_noise": noise,
            "pdf_matrix": pdfMatrix,
            "pdf_bins": pdfBins,
            "max": colMax,
            "min": colMin,
            "mean": colMean,
            "median": colMedian,
            "mode": colMode,
            "nlnm": nlnm,
            "nhnm": nhnm,
            "percent_above_nhnm": pctAbove,
            "percent_below_nlnm": pctBelow,
            "percent_95": col95,
            "percent_90": col90,
            "percent_10": col10,
            "percent_5": col05,
        }
        data.append(d)

    # Insert data into database if available
    if database is not None:
        database.insert_metric(data)

    return data
