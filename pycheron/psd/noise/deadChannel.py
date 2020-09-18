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

__all__ = ["DDT", "isDC"]

import numpy as np
import obspy


def DDT(tr, detrend=True, demean=True, taper=0.1):
    """
    Demeans, detrends, and applies a cosine taper to a trace object

    :param tr: obspy trace object
    :type tr: obspy.core.trace.Trace
    :param detrend: boolean specifying whether to detrend
    :type detrend: bool
    :param demean: boolean specifying whether to demean
    :type demean: bool
    :param taper: proportion (percentage) of the signal to be tapered at each end
    :type taper: float

    :return: trace object that has been demeaned, detrended, and tapered
    :rtype: obspy.core.trace.Trace

    **Example**

    .. code-block:: python

        #import function
        from psd.noise.deadChannel import DDT

        # Initialize IRIS client
        client = Client("IRIS")

        #Grab out P-wave onset for a big earthquake
        starttime = UTCDateTime("2010-02-27T06:30:00.000")
        endtime = UTCDateTime("2010-02-27T07:30:00.000")

        # Grab data from IRIS web client server with specified starttime/endtime and SNCL
        st = client.get_waveforms("IU","ANMO","00","BHZ", starttime,endtime)

        # Get trace from stream
        tr = st[0]

        # Copy stream so doesn't permanently change trace object
        tr1 = tr.copy()

        # Detrend, demean, taper trace
        trDDT = DDT(tr1, True, True, 0.1)

    **Plotting**

    .. code-block:: python

        # Compare detrended, demeaned, and tapered trace to original data
        fig = plt.figure(figsize=(10,10))
        fig.suptitle('Demean-Detrend-Cosine Taper')

        ax1 = fig.add_subplot(311)
        ax1.plot(tr.times(), tr.data, label='Raw Data')
        ax1.hlines(y=0,xmin=np.min(tr.times()), xmax=4000, colors='r', linestyles='solid', label='Raw Data Baseline y= 0')
        ax1.legend(loc = 2)

        ax2 = fig.add_subplot(312)
        ax2.plot(tr.times(), trDDT, label = 'DDT Data')
        ax2.hlines(y=0,xmin=np.min(tr.times()), xmax=4000, colors='r', linestyles='solid', label='DDT Data Baseline y = 0')
        ax2.legend(loc = 2)

        ax3 = fig.add_subplot(313)
        ax3.plot(tr.times(), tr.data-trDDT, label = 'Raw-DDT')
        ax3.hlines(y=0,xmin=np.min(tr.times()), xmax=4000, colors='r', linestyles='solid', label='Raw-DDT Data Baseline y = 0')
        ax3.legend(loc = 2)

    .. image:: _static/ddtPlot.png


    """

    # Handles merged traces
    while True:
        try:
            # demean and detrend
            if demean and not detrend:
                tr.detrend("demean")

            elif detrend and demean:
                # tr.detrend('linear')
                pass

            elif detrend and not demean:
                tr.detrend("simple")

            else:
                tr = tr

            # add taper
            if taper > 0:
                tr.taper(type="cosine", max_percentage=taper)

            if isinstance(tr, obspy.core.trace.Trace):
                return tr
            elif isinstance(tr, obspy.core.stream.Stream) and len(tr) == 1:
                trS = tr[0]
                return trS

        # If tr has gaps and has been merged
        # TODO: validate this section is implemented correctly with splitting traces and merging back together as
        # TODO: not in original R code
        except NotImplementedError:

            # split trace
            st = tr.split()

            # looping through split traces
            for i in range(len(st)):
                trI = st[i]

                # demean and detrend
                if demean and not detrend:
                    trI.detrend("demean")

                elif detrend and demean:
                    trI.detrend("linear")

                elif detrend and not demean:
                    trI.detrend("simple")
                else:
                    trI = trI

                # add taper
                if taper > 0:
                    trI.taper(type="hann", max_percentage=taper)
            # merging dmeaned, detrended, and tapered traces back together
            trM = st.merge()

            if isinstance(trM, obspy.core.trace.Trace):
                return trM
            elif isinstance(trM, obspy.core.stream.Stream) and len(trM) == 1:
                trS = trM[0]
                return trS


def isDC(tr, detrend=False, demean=True, taper=0):
    """

    Tests trace to see whether it is a dead channel, e.g., flatlined (std is zero)

    :param tr: trace object
    :type tr: obspy.core.trace.Trace
    :param detrend: boolean specifying whether to detrend
    :type detrend: bool
    :param demean: boolean specifying whether to demean
    :type demean: bool
    :param taper: proportion (percentage) of the signal to be tapered at each end
    :type taper: float

    :return: True or False as to whether the trace is Dead Channel.
    :rtype: bool

    **Example**

    .. code-block:: python

        #import function
        from pycheron.psd.noise.deadChannel import isDC

        #test data with multiple traces
        data = 'test/test_data/6e_sp06_Eall.397679.tar.mseed'

        #reading in stream
        st = obspy.read(data)

        #grabbing trace
        tr = st[0]

        #test if dead channel
        isDC(tr)
        >>> False

    """
    # demean, detrend, and taper data
    tr = DDT(tr, detrend, demean, taper)
    # initialize rd
    rd = 0
    # rounding to 4 decimal places
    if isinstance(tr, obspy.core.trace.Trace):
        rd = np.round(tr.data, 4)
    if isinstance(tr, obspy.core.stream.Stream):
        try:
            rd = np.round(tr[0].data, 4)
        except IndexError:
            return True

    if np.count_nonzero(rd) == 0:
        return True
    else:
        return False
