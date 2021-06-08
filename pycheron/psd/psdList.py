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

from concurrent.futures import TimeoutError

__all__ = ["psdList"]

from pycheron.psd.McNamaraPSD import McNamaraPSD
from pycheron.util.logger import Logger
from pycheron.util.getChannel import filterChannel
import numpy as np
import sys

from pebble import ProcessPool


def psdList(st, processes=7, logger=None):
    """

    Creates a list of Power Spectral Densities (PSDs) from an ObsPy Stream object.
    Uses the McNamara method [.[#]] to create smoothed (binned) PSDs (uncorrected). The McNamara method breaks up data
    into chunks depending on the channel identifiers with 50% overlap
    (https://cran.r-project.org/web/packages/IRISSeismic/IRISSeismic.pdf).

    :param st: Obspy Stream
    :type st: obspy.core.stream.Stream
    :param processes: Number of processes to create. **Note:** It is advised you do not use more than 7 processors if
        on a personal computer.
    :type processes: int
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger

    :return: List of psds, with the following access flags:

             * ``psd[i][j][0]`` = frequency
             * ``psd[i][j][1]`` = psd
             * ``psd[i][j][2]`` = snclq
             * ``psd[i][j][3]`` = Starttime of data segment
             * ``psd[i][j][4]`` = Starttime of data segment

             where ``i`` = which trace in stream, ``j`` = which data segment

    :rtype: list

    * Code originally ported from IRISSeismic R Cran Package
      (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
      The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index) and
      augmented and adapted for use within Pycheron

    **Example**

    .. code-block:: python

        # import function
        from pycheron.psd.psdList import psdList

        # test data
        data1 = 'test/test_data/7a_cabn_bhe.884965.tar.mseed'

        # reading in stream
        st = obspy.read(data1)

        # print out stream information
        print st

        # creating list of psds (uncorrected)
        psd = psdList(st)

        # The output is a nested list. The first list (psd[i]) is a list of all of the different traces in the stream.
        # Since the sample data only has one trace, the list has a length of one. The second nested list (psd[i][j])
        # is a list of all of the shorter segments, which the function broke it up into. The last list
        # (psd[i][j][k]) is the individual slice of the  data with frequency, psd, snclq, starttime and endtime.

        # Output not show for brevity
        print 'Frequency:', psd[0][0][0] #Frequency of first segment or first trace
        print 'PSD:', psd[0][0][1] #PSD of first segment or first trace
        print 'SNCLQ:', psd[0][0][2] #SNCLQ of first segment or first trace
        print 'Starttime:', psd[0][0][3] #Starttime of data segment
        print 'Endtime:', psd[0][0][4] #Endtime of data segment

    .. rubric:: References

    .. [#] McNamara, D. E., & Boaz, R. I. (2006). Seismic noise analysis system using power spectral density probability
           density functions: A stand-alone software package. US Geological Survey.

    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Filters out SOH channel traces from stream
    stN = filterChannel(st)
    # If gaps, merge traces -- first create copy of stream so doesn't permanently change stream object
    if len(stN.get_gaps()) != 0:
        stM = stN.copy()
        # Loop through traces in stream object and round sampling rates -- this is to get around the fact that for small
        # deviations in sampling rate, e.g., 50 and 49.9999 ObsPy won't merge the traces. For large deviations this will
        # throw an error on merge, but for small ones this should be fine
        for i in range(len(stM)):
            tr1 = stM[i]
            tr1.stats.sampling_rate = np.round(tr1.stats.sampling_rate, 2)
        stM = stM.merge(fill_value=0)
    else:
        stM = stN.copy()

    # Initialize output list
    channel = []
    PSDList = []
    # If length greater than 1, start a pool
    if len(stM) > 1:
        with ProcessPool(max_workers=processes) as pool:
            p = pool.map(_pool_psd_python, stM, timeout=600)
            r = p.result()
            while True:
                try:
                    PSDList.append(next(r))
                except StopIteration:
                    break
                # If timeout, exit
                except TimeoutError as error:
                    logger.error("psdList(): Python psdList Timed Out" + error.args[0])
                    sys.exit()
            # Let users know calculation complete
            logger.log("psdList(): Python PSD Complete")

        # Loop through stream and append tr.stats.channel to channel list
        for i in range(len(stM)):
            tr = stM[i]
            channel.append(tr.stats.channel)
    # If it's not, initalize list, start pool and append to list
    else:
        PSDList = []
        psd = _pool_psd_python(stM[0])
        tr = stM[0]
        PSDList.append(psd)
        channel.append(tr.stats.channel)

    return PSDList


def _pool_psd_python(tr):
    # If data not at least an hour return
    if tr.stats.endtime - tr.stats.starttime < 3600:
        return

    # Create temp to fill psds into output for every run
    temp = []
    channel = tr.stats.channel
    # Set high frequency boundary, which will always be set to Nyquist no matter the channel band code is
    hiFreq = 0.5 * tr.stats.sampling_rate

    # Set center frequency for octaves to be aligned, e.g., alignment frequency from which octaves will be generated
    alignFreq = 0.1

    # "C"hoose chunk size based on channel band code (=sampling rate). See:
    # http://www.fdsn.org/seed_manual/SEEDManual_V2.4.pdf , Appendix A"
    # (Callahan, 2020)
    if channel.startswith("L"):
        z = 3 * 3600
        loFreq = 0.001

    # V
    elif channel.startswith("V"):
        z = 24 * 3600
        loFreq = 0.0001
        # Special for VH
        alignFreq = 0.025
    # M
    elif channel.startswith("M"):
        z = 2 * 3600
        loFreq = 0.0025

    # Changed to z = 10, loFreq = 0.6 from IRISSeismic 1.4.5 to 1.4.6, but keeping at 1.4.5 options for now because
    # 10s increments take too long for processing but do increase the statistical significance -- think about in
    # future
    else:
        z = 3600
        loFreq = 0.005
    # else:
    #    z = 10
    #    loFreq = 0.6

    # initialize a bunch of values, like start, end time, secs, secsLeft
    endtime = tr.stats.endtime
    starttime = tr.stats.starttime
    index = 0
    start = starttime
    end = start + z
    secs = endtime - start
    secsLeft = secs

    # While 99% of the chunk still remains loop through
    while secsLeft >= 0.99 * z:
        # slice trace
        sl = tr.slice(start, end)

        # check if Dead channel
        # if not isDC(sl): //TODO Change this in future to run dead channels becaause won't be calculated anyway
        # calculate PSDs using McNamaraPSD
        psds = McNamaraPSD(sl, loFreq, hiFreq, alignFreq, binned=True)
        temp.append(psds)

        # Increment the window with 50% overlap before continuing on through the while statement
        index += 1
        start += z / 2
        end += z / 2
        secsLeft = tr.stats.endtime - start

    return temp
