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

__all__ = ["triggerOnset"]
import numpy as np
from scipy.stats.mstats import mquantiles


def triggerOnset(x, picker, threshold=None, index=False):
    """

    "Uses data returned by STALTA "first break picking" method and a user selected threshold to determine the arrival
    time of a seismic event." (Callahan, 2020)

    :param x:  Obspy Trace object
    :type x: obspy.core.trace.Trace
    :param picker: Results from applying the STALTA method to this trace
    :type picker: numpy.array
    :param threshold: Optional threshold at which triggering should occur.
                      "Note the appropriate value for the threshold will depend upon the exact STA/LTA algorithm used
                      and the noise level in the signal." (Callahan, 2020)
    :type threshold: float
    :param index: Optional boolean to return the index (rather than the time) of event onset.
    :type index: bool

    :return: Returns a single value identifying the onset of the seismic event or None if nothing is detected. The
             returned value will be a time by default, or a numeric index if index = True.
    :rtype: float or None

    Method identifies the point at which the picker first rises above the threshold, e.g., only returns the first
    time point or None at which STALTA is above threshold. It does not implement the full functionality found in Obspy
    method of the same name as is not currently needed at this point, but maybe a TODO later:
    (http://docs.obspy.org/packages/autogen/obspy.signal.trigger.triggerOnset.html)

    When no threshold is supplied, an appropriate value is calculated from the picker with:

    ``threshold = scipy.stats.mstats.quantile(picker, 0.999)``

    * Code originally ported from IRISSeismic R Cran Package
      (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
      The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index) and
      augmented and adapted for use within Pycheron.

    **Example**

    .. code-block:: python

        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        from pycheron.sigpro.STALTA import STALTA
        from pycheron.sigpro.triggerOnset import triggerOnset

        # Instantiate client object
        client = Client("IRIS")

        # Grab data from 2010-02-27 06:00:00 to 2010-02-27 09:00:00
        t = UTCDateTime("2010-02-27T06:00:00.000")
        st = client.get_waveforms("IU","ANMO","00","BHZ",t,t+180*60)
        tr = st[0]
        picker = STALTA(tr,3,30)

    **Plotting**

    .. code-block:: python

        import matplotlib.pyplot as plt

        # Using results from above example
        plt.plot(tr)
        plt.axvline(x=to,color = 'red')

    .. image:: _static/triggerOnset.png

    """

    # If threshold and index set to None, set appropriately
    # Set threshold
    if threshold is None:
        threshold = mquantiles(picker[~np.isnan(picker)], 0.99999)
    # Set index
    if index is None:
        index = False

    # Find first index above the threshold provided. In case only one, grab only the first index
    eventIndex = int(np.argwhere(picker[~np.isnan(picker)] > threshold)[0])

    # Add indices that are nans into the calculation as otherwise will be too early. Had to ignore them in above
    # calculation as otherwise would just give nan value
    np_end = int(np.argwhere(~np.isnan(picker))[0])
    eventIndex += np_end

    # Return the index where rises above the threshold, if desired
    if index is True:
        return eventIndex

    # Otherwise calculate the time at which that index occurs and return it instead
    if not eventIndex:
        eventTime = None
    else:
        eventTime = x.stats.starttime + eventIndex / x.stats.sampling_rate
    return eventTime
