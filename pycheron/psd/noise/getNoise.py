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

__all__ = ["getNoise"]

from obspy.clients.iris.client import Client
from obspy import UTCDateTime
from ispaq.evalresp import evalresp as read_Evalresp
from pycheron.util.logger import Logger
from pycheron.dataAcq.css import platform_file_path
from pycheron.dataAcq.fap import fap_reader
import numpy as np


def getNoise(psd_list, evalresp=None, logger=None, units='acc'):
    """
    Obtains instrument correction information at the specified frequencies from IRIS/user provided evalresps.
    The instrument correction is applied to every PSD in the list and the now corrected PSD values are returned
    as a list. This is Pycheron's version of IRISSeismic's psdList2NoiseMatrix function.

    :param psd_list: list of PSDs (uncorrected), use `pycheron.psdList.psdList` function
    :type psd_list: list
    :param evalresp: dict containing tuple of respfile information (file location, file, file type).
                     Either obtained from client.evalresp if None, or can be user supplied.
    :type evalresp: dict {channel (str): resp file information (file location, file, file type)}
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger
    :param units: Units to return response in. Can be either def, dis, vel or acc. If none, specified, will use default units. 
    :type units: string


    :return:

     * List of frequencies for associated psd
     * List of noise matrices for each corrected PSD (for each trace of stream)
     * Modified input list (psd_list), will be the same length if no stations are popped,
       otherwise may be a different length.

    :rtype:

        * `list`
        * `list`
        * `list`

     * Code originally ported from IRISSeismic R Cran Package
       (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
       The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index)
       and augmented and adapted for use within Pycheron

    **Example**

    .. code-block:: python

        # import psdList since psdList is input to getNoise
        from pycheron.psd.psdList import psdList
        from pycheron.psd.noise.getNoise import getNoise

        # test data
        data = 'test/test_data/7a_cabn_bhe.884965.tar.mseed'

        # reading in stream
        st = obspy.read(data)

        # getting noise to input into function
        # creating list of psds, this is a REQUIREMENT as the psdStatistic takes in a psdList
        psds = psdList(st)
        # Get instrument corrected psds
        f,n,psd = getNoise(psds)

    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Initialize a few variables
    client = Client(timeout=30)
    noiseData = []
    freqData = []
    indices = []

    # Ensure units are properly supplied, otherwise log an error
    if units in ('acc', 'def', 'dis', 'vel'):
        pass
    else:
        logger.error("getNoise(): units are not defined properly for response removal. Valid units are 'ACC', 'DEF', 'DIS', 'VEL'")

    # First loops through psd_list and removes empty psds
    pop_i = []
    for i in range(len(psd_list)):
        # If psd is empty, it will throw error
        try:
            len(psd_list[i])
        except TypeError:
            # adds empty indices to list
            pop_i.append(i)

    # Pop indices
    psd_list = _pop_indices(pop_i, psd_list)

    # Loop through PSDs in psd_list
    for i in range(len(psd_list)):
        # Get a count
        psdCount = len(psd_list[i])
        # parsing psd_list for individual PSDs
        # Obtain snclq information needed for evalresp from psd_lists
        try:
            network, station, location, channel, quality = str(psd_list[i][0][2][0]).split(".")
            starttime = psd_list[i][0][3]
            freq = psd_list[i][0][0]
        except ValueError:
            network, station, location, channel = str(psd_list[i][0][2][0]).split(".")
            starttime = psd_list[i][0][3]
            freq = psd_list[i][0][0]

        # Get min, max, and number of freqs
        minFreq = np.nanmin(freq)
        maxFreq = np.nanmax(freq)
        nFreq = len(freq)

        # Get instrument response from IRIS webservice in the case where evalresp = None
        if evalresp is None:
            try:
                iris = client.evalresp(
                    network,
                    station,
                    location,
                    channel,
                    starttime,
                    minFreq,
                    maxFreq,
                    nFreq,
                    units,
                    output="fap",
                )
                amp = []
                for j in range(len(iris)):
                    amp.append(iris[j][1])
            except Exception as ex:
                logger.error("getNoise(): " + str(ex))
                indices.append(i)
                logger.error(
                    "getNoise(): Popping off element {0}.{1}.{2}.{3}".format(
                        str(network), str(station), str(location), str(channel)
                    )
                )
                continue
        # If evalresp not none, use user specified evalresp input, either FAP or EVRESP
        else:
            e_filtered = evalresp[channel]
            respfile = platform_file_path(filename=e_filtered[1], dirname=e_filtered[0])

            if e_filtered[2].upper() == "FAP":
                iris = fap_reader(respfile)

            elif e_filtered[2].upper() == "EVRESP":
                try:
                    iris = read_Evalresp(
                        sfft=minFreq,
                        efft=maxFreq,
                        nfft=nFreq,
                        filename=respfile,
                        date=UTCDateTime(starttime),
                        station=station,
                        channel=channel,
                        network=network,
                        locid=location,
                        units=units,
                        debug=True,
                    )
                    amp = iris[1]
                except Exception as e:
                    raise ValueError(f"Failed to read EVRESP file in getNoise: {e}")

            amp = iris[1]

        # Initialize amplitude container, then append to it

        # Convert instrument response to dB and square it
        correction = 10 * np.log10(amp) * 2

        # Create empty rawNoiseMatrix to fill based on length of psd_list and number of frequencies
        rawNoiseMatrix = np.empty([psdCount, nFreq])

        # Loop through list of psds and fill raw noise matrix
        for k in range(psdCount):
            try:
                rawNoiseMatrix[k] = psd_list[i][k][1]
            # Set index to nan if hit a value error
            except ValueError:
                rawNoiseMatrix[k] = np.nan

        # Construct correction Matrix that matches dimensions of noiseMatrix
        correctionMatrix = np.tile(correction, psdCount)
        correctionMatrix = np.reshape(correctionMatrix, (psdCount, nFreq))

        # Divisiion in dB space is just subtraction. We want to divide noise matrix by correction matric, so we subtract
        noiseMatrix = rawNoiseMatrix - correctionMatrix

        # Append to appropriate containers then return
        noiseData.append(noiseMatrix)
        freqData.append(freq)

    # pop out indices that are bad (e.g., returned HTTP errors above)
    psd_list = _pop_indices(indices, psd_list)

    return freqData, noiseData, psd_list


def _pop_indices(indices, arr):
    """
    Pops values from array by indices
    :param indices: List of indices to pop
    :param arr: Array containing values that need to be removed
    :return: Array with removed values
    """

    # Initialize count
    count = 0
    # Loop through provided indices
    for i in range(len(indices)):
        if i == 0:
            arr.pop(indices[i])
            count += 1
        else:
            # calculating new position of index since the shifted.
            new_i = indices[i] - count
            arr.pop(new_i)
            count += 1
    return arr
