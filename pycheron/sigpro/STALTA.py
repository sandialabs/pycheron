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

__all__ = ["STALTA"]

from pycheron.psd.noise.deadChannel import DDT
from pycheron.rollseis.roll_mean import roll_mean
from pycheron.rollseis.roll_stalta import roll_stalta
from pycheron.sigpro.envelope import envelope
from pycheron.util.logger import Logger


def STALTA(
    tr,
    staSecs=3,
    ltaSecs=30,
    algorithm="classic_LR",
    increment=1,
    demean=True,
    detrend=True,
    taper=0,
    logger=None,
    fortran=False,
):
    """

    "The STALTA method applies one of several STA/LTA "first break picking" algorithms (``classic LR``, ``classic RR``,
    or ``EarleAndShearer``) to ObsPy Trace data to automatically detect seismic events." (Callahan, 2020)

    :param tr: ObsPy Trace object
    :type tr: obspy.core.trace.Trace
    :param staSecs: length of the short averaging window in secs (default = 3)
    :type  staSecs: int
    :param ltaSecs: length of long averaging window in secs (default = 30)
    :type ltaSecs: int
    :param algorithm: STA/LTA algorithm to use. Options include ``classic_LR``, ``classic_RR``, and
                      ``EarleAndShearer_envelope`` (default = ``classic_LR``)
    :type algorithm: str
    :param increment: the increment to use when sliding the averaging windows to the next location.
    :type increment: str
    :param demean: boolean flag determining whether to demean the data before applying the algorithm
    :type demean: bool
    :param detrend: boolean flag determining whether to detrend the data before applying the algorithm
    :type detrend: bool
    :param taper: proportion of the signal to be tapered at each end before applying the algorithm (default = 0)
    :type taper: int
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger
    :param fortran: Use Fortran libs or not. If libs will not compile or on a Windows Machine, set to False
    :type fortran: bool

    :return: Matrix of values of the same length as tr.data. "Note the returned vector will contain
             NaNs near the edges of the trace where insufficient data are available to fill windows. Additional NaN
             values will appear for every index that is skipped over when the increment parameter is greater than 1."
             (Callahan, 2020)
    :rtype: numpy.array


    "By default, this method uses the ``"classic_LR"`` algorithm which calculates the average power in the trace data
    over a short window (STA) and a long window (LTA). With this algorithm, windows are "left/right aligned" meaning
    that the point for which STA/LTA is calculated is at the lefttmost edge of the STA window and the rightmost edge
    of the LTA window. The resulting STA/LTA ratio thus has the same number of points as the original data. This is a
    standard method of "first break picking" and can be used to identify the onset of a seismic event. [#]_ [#]_ [#]_

    **Methodology**

    1. ``algorithm = "classic_RR"``. This is the original STA/LTA algorithm with "right alignment" [#]_

    .. math::

        STA(x_{i}) = (1/ns)\sum_{j = i-ns}^{i} x_{i}^{2}

        LTA(x_{i}) = (1/nl)\sum_{j = i-nl}^{i} x_{i}^{2}

        r_{i} = STA_{i}/LTA_{i}

    .. code-block:: console

       [---------- LTA ---------*]
                       [-- STA -*]

    2. ``algorithm = "classic_LR"`` (default). This algorithm has the index at the left edge of the STA window and the
       right edge of the LTA window. [#]_

    .. math::

        STA(x_{i}) = (1/ns)\sum_{j = i+ns}^{i} x_{i}^{2}

        LTA(x_{i}) = (1/nl)\sum_{j = i-nl}^{i} x_{i}^{2}

        r_{i} = STA_{i}/LTA_{i}

    .. code-block:: console

       [---------- LTA --------*]
                              [*- STA --]

    3. ``algorithm = "EarleAndShearer_envelope"`` [#]_

    .. math::

       STA(x_{i}) = (1/ns)\sum_{j=i}^{i+ns} mod(H(x))_{i}

       LTA(x_{i}) = (1/nl)\sum_{j=i-nl}^{i} mod(H(x))_{i}

       r_{i} = STA_{i}/LTA_{i}

    where: *H(x)* is the hilbert transform of the data and *mod(H(x))* is the 'envelope' of the seismic signal

    .. code-block:: console

       [---------- LTA ---------*]
                               [*- STA --]" (Callahan, 2020)

    .. note:: The EarleAndShearer algorithms uses Hilbert transforms which involves performing an FFT, thus it can
              take significantly longer than the "classic" algorithms for longer seismic signals (>500,000 pts)}
              (Callahan, 2020). "For higher resolution channels, picking an increment of 2/sampling rate can greatly
              speed up processing times and still generate reasonable results." (Callahan, 2020)

    * Code originally ported from IRISSeismic R Cran Package
      (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
      The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index) and
      augmented and adapted for use within Pycheron.

    **Example**

    .. code-block:: python

        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        from pycheron.sigpro.STALTA import STALTA

        #Instantiate client object
        client = Client("IRIS")

        #Grab data from 2010-02-27 to 2010-02-28
        t = UTCDateTime("2010-02-27T00:00:00.000")
        st = client.get_waveforms("IU","ANMO","00","BHZ",t,t+1440*60)
        tr = st[0]
        picker = STALTA(tr,3,30)

    **Plotting**

    .. code-block:: python

        import matplotlib.pyplot as plt

        # Plot the trace and overlay the picker (on separate axes for easier viewing)
        fig, ax1 = plt.subplots()
        # Uses output from above example
        ln1 = ax1.plot(tr,color='blue', label = 'raw data')
        ax1.set_xlabel('Number of points')
        ax1.set_ylabel('Amplitude', color ='blue')

        # Create picker axis and plot
        ax2=ax1.twinx()
        # Uses output from above example
        ln2 = ax2.plot(picker, color = 'red', label = 'picker')
        ax2.set_ylabel('Amplitude', color = 'red')

        # Create combined legend
        lns = ln1+ln2
        labs = [l.get_label() for l in lns]
        ax1.legend(lns,labs)

    .. image:: _static/stalta.png

    .. rubric:: References

    .. [#] First break picking: multi-purpose STA/LTA trigger algorithm, http://en.wikipedia.org/wiki/First_break_picking
    .. [#] "A Comparison of Select Trigger Algorithms for Automated Global Seismic Phase and Event Detection"
            http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.116.245&rep=rep1&type=pdf
    .. [#] Wong, "Automatic time-picking of first arrivals on noisy microseismic data"
           http://www.crewes.org/ForOurSponsors/ConferenceAbstracts/2009/CSEG/Wong_CSEG_2009.pdf
    .. [#] Wong, GEOPHYSICS, VOL. 75, NO. 4 (JULY-AUGUST 2010); P. V67-V76 Automatic first-breaks picking: New
           strategies and algorithms". http://www.fcaglp.unlp.edu.ar/~velis/papers/PickingGeop10.pdf
    .. [#] Wong, "Adaptive microseismic event detection and automatic time picking"
           http://www.cspg.org/documents/Conventions/Archives/Annual/2012/279_GC2012_Adaptive_Microseismic_Event_Detection.pdf
    .. [#] Earle and Shearer, "Characterization of Global Seismograms Using an Automatic-Picking Algorithm" Bulletin of the
           Seismological Society of America, Vol. 84, No. 2, pp. 366-376, April 1994


    """
    # TODO: We should probably handle errors where ltaSecs or increment are
    # passed in as 0 (or negative values), currently ltaSecs will cause a divide
    # by zero error if passed as zero. staSecs has some error handling below.

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # Calculate STA window and LTA windows given the length of the window and sampling rate of the data
    n_sta = int(staSecs * tr.stats.sampling_rate)
    n_lta = int(ltaSecs * tr.stats.sampling_rate)

    # Test for windows that are too short based on sampling rate
    if n_sta < 1:
        logger.error(
            "STALTA(): STA window of %s secs is too short for a trace sampling rate of %s Hz"
            % (staSecs, tr.stats.sampling_rate)
        )
        return

    # If demean/detrend True, demean, detrend data, make copy of trace before though so we don't permanently change
    # original trace
    if demean or detrend:
        trC = tr.copy()
        trC = DDT(trC, demean, detrend, taper)
    # Otherwise trC is just the original data
    else:
        trC = tr

    # Check that there is sufficient data
    if len(trC) < n_lta:
        logger.error("STALTA(): insufficient data")
        return

    # DC signals will have all zeroes after demeaning. Return all zeroes in this case. May want to use isDC here?
    if trC.data.all() == 0:
        return trC.data

    # Apply the user specific algorithm to the data
    # Classic RR algorithm
    if algorithm == "classic_RR":
        # For classic RR, LTA/STA windows use the same right edge. Right-Right means that the index is in both windows
        if increment != 1:
            logger.error("STALTA(): algorithm classic_RR requires increment = 1")
            return
        # Calculate rolling mean STA, LTA, then get first breaks
        sta = roll_mean(trC.data ** 2, n_sta, increment, align="right")
        lta = roll_mean(trC.data ** 2, n_lta, increment, align="right")
        picker = sta / lta
    # This is the same algorithm as the classic RR, but with left-right alignment. Left-right alignment means the index
    # is in both windows
    # Classic LR algorithm
    elif algorithm == "classic_LR":
        # Get first breaks
        picker = roll_stalta(trC.data ** 2, n_sta, n_lta, increment, fortran, logger)
    # Earle and Shearer algorithm
    elif algorithm == "EarleAndShearer_envelope":
        # Error out if increment isn't 1
        if increment != 1:
            logger.error(
                "STALTA(): algorithm EarleAndShearer_envelope requires increment = 1"
            )
            return
        # Get the enveloper
        data = envelope(trC)
        # Calculate the rolling mean
        sta = roll_mean(data, n_sta, increment, align="left")
        lta = roll_mean(data, n_lta, increment, align="right")
        # Get the first breaks
        picker = sta / lta
    # Throw error if user specified algorithm that isn't available
    else:
        logger.error("STALTA(): algorithm %s not recognized" % algorithm)
        return

    return picker
