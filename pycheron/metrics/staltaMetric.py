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

__all__ = ["staltaMetric"]

import numpy as np
from pycheron.sigpro.STALTA import STALTA
import multiprocessing as mp
import time
from pycheron.util.logger import Logger


def staltaMetric(
    st,
    staSecs=3,
    ltaSecs=30,
    increment=1,
    algorithm="classic_LR",
    logger=None,
    processes=1,
    fortran=True,
    database=None,
):
    """
    Calculates the maximum STA/LTA for a stream of seismic data

    :param st: obspy stream object
    :type st: obspy.core.stream.Stream
    :param staSecs: length of short term averaging window in secs
    :type staSecs: int
    :param ltaSecs: length of long term averaging window in secs
    :type ltaSecs: int
    :param increment: increment (in secs) used when sliding the average windows to the next location
    :type increment: int
    :param algorithm: algorithm to be used. Currently supported algorithms include: "classic_RR", "classic_LR",
                        "EarleAndShearer_envelope"
    :type algorithm: str
    :param logger: If using a logger, (you must create one using the util.logger class)
    :type logger: pycheron.util.logger.Logger
    :param processes:  Number of processes to use for calculation (DEFAULT = 3)
    :type processes: int
    :param fortran: Use Fortran libs or not. If libs will not compile or on a Windows Machine, set to False
    :type fortran: bool
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database

    :return: list of dictionaries containing the following keys and types:

                * snclq (`str`)
                * start_time (`str`)
                * end_time (`str`)
                * max_stalta (`float`)
                * event_time (`str`)
                * metric_name (`str`)

    :rtype: list

    **Algorithm Notes**

    This metric applies the STALTA method to every trace in the stream object using the following parameter settings:

    * demean = True
    * detrend = True
    * taper = 0.0

    The final metric value is the maximum STALTA value found in any trace in the stream object.

    **Example**

    .. code-block:: python

        # Initialize IRIS client
        client = Client("IRIS")

        # # Define start/end time and then get stream object
        starttime = UTCDateTime("2012-12-12T00:00:00.000")
        endtime = UTCDateTime("2012-12-13T00:00:00.000")
        st = client.get_waveforms("AK","GHO","","BHN",starttime, endtime)

        #Calculate the STA/LTA metric (use all default parameters) and show the results

        stalta = staltaMetric(st)
        print stalta
        >>> [{'metric_name': 'staltaMetric', 'event_time': "2012-12-12T21:37:25.168400Z", 'start_time': "2012-12-12T00:00:00.008400Z" , 'snclq': u'AK.GHO..BHNM', 'max_stalta': 75.67122, 'end_time': "2012-12-12T23:59:59.988400Z"}]

    .. rubric:: References

    .. [#] http://en.wikipedia.org/wiki/First_break_picking (Wikipedia)
    .. [#] http://www.crewes.org/ForOurSponsors/ConferenceAbstracts/2009/CSEG/Wong_CSEG_2009.pdf (Wong et. al. 2009)
    .. [#] http://www.fcaglp.unlp.edu.ar/~velis/papers/PickingGeop10.pdf (Sabbione and Velis 2010)

    """
    if logger == None:
        logger = Logger(None)

    # If gaps, merge, but first make a copy of the data so don't overwrite original contents (either by merging gaps
    # and/or demeaning/detrending within STALTA)
    st1 = st.copy()

    if len(st1.get_gaps()) != 0:
        logger.warn("staltaMetric(): This data has gaps; merging traces...")
        # Copy stream so don't overwrite completely, then merge traces with gaps
        st1 = st1.merge()

    args = []
    for i in range(len(st1)):
        args.append([st1[i], staSecs, ltaSecs, increment, algorithm, logger, fortran])
    if len(st1) > 1:
        pool = mp.Pool(processes=processes)
        out = pool.map(_stalta_pool, args)
        pool.close()

    else:
        out = []
        d = _stalta_pool(args[0])
        out.append(d)

    if database is not None:
        database.insert_metric(out)

    return out


def _stalta_pool(args):
    start = time.time()
    maxSTALTA = 0
    tr = args[0]
    staSecs = args[1]
    ltaSecs = args[2]
    increment = args[3]
    algorithm = args[4]
    logger = args[5]
    fortran = args[6]

    starttime = tr.stats.starttime
    endtime = tr.stats.endtime

    # Initialize eventTime
    eventTime = starttime

    # Testing for enough data
    nlta = ltaSecs * tr.stats.sampling_rate
    nsta = staSecs * tr.stats.sampling_rate
    if len(tr) <= (nlta + nsta):
        logger.error(
            "staltaMetric(): Not enough data for trace %s" % (str(tr.get_id()))
        )
        return

    stalta = STALTA(
        tr,
        staSecs,
        ltaSecs,
        algorithm,
        increment,
        demean=True,
        detrend=True,
        taper=0,
        logger=logger,
        fortran=fortran,
    )
    # Add check for when result blows up to infinity, replace with NaNs

    (rinfs,) = np.where(np.isinf(stalta) is True)
    for j in range(len(rinfs)):
        stalta[rinfs[j]] = np.nan

    # Check for scenario where entire result is NaNs
    if np.isnan(stalta).all():
        logger.error(
            "staltaMetric(): stalta returns a vector with all NaNs for trace %s"
            % (str(tr.get_id()))
        )
        return

    # Find Maxiumum STA/LTA value
    traceMaxSTALTA = np.nanmax(stalta)

    # Calculate time at which STA/LTA max occurred
    (eventIndex,) = np.where(stalta == traceMaxSTALTA)
    traceEventTime = starttime + eventIndex[0] / tr.stats.sampling_rate

    if traceMaxSTALTA > maxSTALTA:
        maxSTALTA = traceMaxSTALTA
        eventTime = traceEventTime
    # if tr is not from mseed
    try:
        snclq = tr.get_id() + "." + tr.stats.mseed.dataquality
    except AttributeError:
        snclq = tr.get_id()

    d = {
        "snclq": snclq,
        "start_time": starttime.isoformat(),
        "end_time": endtime.isoformat(),
        "max_stalta": maxSTALTA,
        "event_time": eventTime.isoformat(),
        "metric_name": "staltaMetric",
    }
    timestart = (time.time() - start) / 60
    logger.log("staltaMetric(): Time in minutes: " + str(timestart))
    return d
