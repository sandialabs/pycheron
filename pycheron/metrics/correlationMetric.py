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

__all__ = ["correlationMetric"]

import numpy as np
from pycheron.psd.noise.deadChannel import isDC
from pycheron.util.logger import Logger
from scipy.stats.stats import pearsonr


def correlationMetric(tr1, tr2, logger=None, database=None):
    """
    Calculates the correlation between two traces of seismic data from same location.

    :param tr1: obspy trace object 1
    :type tr1: `obspy.core.trace.Trace`
    :param tr2: obspy trace object 2
    :type tr2: `obspy.core.trace.Trace`
    :param logger: If using a logger, (you must create one using the util.logger class)
    :type logger: pycheron.util.logger.Logger
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database

    :return: Dictionary containing the following keys and types:

            * snclq (`str`)
            * start_time (`str`)
            * end_time (`str`)
            * correlation_coefficient (`float`)
            * p_value (`float`)
            * metric_name (`str`)

    :rtype: dict

    The correlation returned is a value in the range [+/-1]. This 'pearson r' correlation is a measure of the strength
    and direction of the linear relationship between two variables that is defined as the (sample) covariance of the
    variables divided by the product of their (sample) standard deviations. Pearson's correlation requires that each
    dataset be normally distributed. Values may vary between -1 and +1 with 0 implying no correlation. Correlations of
    -1 or +1 imply an exact linear relationship. Positive correlations imply that as x increases, so does y. Negative
    correlations imply that as x increases, y decreases.

    The p-value roughly indicates the probability of an uncorrelated system producing datasets that have a Pearson
    correlation at least as extreme as the one computed from these datasets. The null hypothesis is that the two
    variables are uncorrelated. The p-value is a number between zero and one that represents the probability that your
    data would have arisen if the null hypothesis were true.The p-values are not entirely reliable but are probably
    reasonable for datasets larger than 500 or so (majority of above statements taken from scipy manpage)

    Seismic traces passed to `correlationMetric` must have the same network and station, must cover the same time range
    and must have the same sampling rate

    **Examples**

    .. code-block:: python

        # Initialize client object
        client = Client("IRIS")
        # Get a days worth of data from 03/01/2013-03/02/2013 and pull out individual traces for comparison
        t = UTCDateTime("2013-03-01T00:00:00.000")
        stANMO = client.get_waveforms("IU","ANMO", "00","BH*",t,t+1440*60)
        tr1 = stANMO[0]
        tr2 = stANMO[1]
        # Compare traces; BH1:BH2
        corr12 = correlationMetric(tr1,tr2)
        print corr12
        >>> {'p-value': 0.0, 'metric_name': 'correclationMetric', 'start_time': UTCDateTime(2013, 3, 1, 0, 0, 0, 19500), 'correlation_coefficient': -0.2041723553359387, 'snclq': u'IU.ANMO.00.BH1.M:IU.ANMO.00.BH2.M', 'end_time': UTCDateTime(2013, 3, 1, 23, 59, 59, 969500)}

    """

    if logger == None:
        logger = Logger(None)

    # -----------Compatibility Checks------------------------------------------
    # Start time and end time for trace data 1
    starttime = tr1.stats.starttime
    endtime = tr1.stats.endtime

    # test if start times same for both traces
    if tr2.stats.starttime != starttime:
        logger.error("correlationMetric(): Start times do not match.")
        return

    # test if end times same for both traces
    if tr2.stats.endtime != endtime:
        logger.error("correlationMetric(): End times do not match.")
        return

    # test if sampling rate same for both traces
    if tr1.stats.sampling_rate != tr2.stats.sampling_rate:
        logger.error("correlationMetric(): Sampling Rates do not match.")
        return

    # Correlation demands vectors of the same length (even though we will ignore NaN values). We will truncate by up to
    # one second to ensure this.
    # Only complain if the sample lengths differ by more than two samples
    l1 = len(tr1)  # Length of trace data 1
    l2 = len(tr2)  # Length of trace data 2
    if abs(l2 - l1) > 2:
        logger.error(
            "correlationMetric(): Incompatible lengths tr1 = %s, tr2 = %s" % (l1, l2)
        )
        return
    else:
        min_length = min(l1, l2)

    # Test that network and station match:
    if tr1.stats.network != tr2.stats.network or tr1.stats.station != tr2.stats.station:
        logger.error(
            "correlationMetric(): Incompatible trace ids, tr1 = %s, tr2 = %s"
            % (tr1.get_id(), tr2.get_id())
        )
        return

    # Check for flat-lined data
    # Make a copy of tr1, tr2 as isDC applies tapering and demeaning techniques, etc that permanently change the trace
    tr1c = tr1.copy()
    tr2c = tr2.copy()
    if isDC(tr1c) is True:
        logger.error(
            "correlationMetric(): %s has one unique sample value (flat-lined). Standard deviation is "
            "zero, correlation is undefined." % (tr1.get_id())
        )
        return

    if isDC(tr2c) is True:
        logger.error(
            "correlationMetric(): %s has one unique sample value (flat-lined). Standard deviation is "
            "zero correlation is undefined." % (tr2.get_id())
        )
        return

    # Merging snclq trace information, station, network will be the same, as likely will be mseed quality, but decided
    # to merge entire snclq information instead of just merging location and channel information like R code does
    try:
        snclq = (
            tr1.get_id()
            + "."
            + tr1.stats.mseed.dataquality
            + ":"
            + tr2.get_id()
            + "."
            + tr2.stats.mseed.dataquality
        )
    except AttributeError:
        snclq = tr1.get_id() + ":" + tr2.get_id()

    # -----------Correlation Function------------------------------------------
    # Calculate correlation metric and P-value
    # Prior to calculation ignore NaNs (if present)
    tr1 = tr1.data[~np.isnan(tr1.data)]
    tr2 = tr2.data[~np.isnan(tr2.data)]
    corr = pearsonr(tr1[:min_length], tr2[:min_length])
    c = {
        "snclq": snclq,
        "start_time": starttime,
        "end_time": endtime,
        "correlation_coefficient": corr[0],
        "p_value": corr[1],
        "metric_name": "correlationMetric",
    }

    if database is not None:
        database.insert_metric(c)

    return c
