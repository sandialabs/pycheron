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
from pycheron.db.sqllite_db import Database
from scipy.stats.stats import pearsonr


def correlationMetric(tr1, tr2, logger=None, database_config=None):
    """
    Calculates the correlation between two ObsPy Traces of seismic data from same location.

    :param tr1: ObsPy Trace object 1
    :type tr1: `obspy.core.trace.Trace`
    :param tr2: ObsPy Trace object 2
    :type tr2: `obspy.core.trace.Trace`
    :param logger: If using a logger, (you must create one using the util.logger class)
    :type logger: pycheron.util.logger.Logger
    :param database_config: dictionary containing the necessary parameters to create
                            a pycheron Database object. 
                            These include "db_name", "session_name", "overwrite", "manual", "wfdb_conn"
    :type database_config: dict

    :return: Dictionary containing the following keys and types:

            * snclq (`str`)
            * start_time (`str`)
            * end_time (`str`)
            * correlation_coefficient (`float`)
            * p_value (`float`)
            * metric_name (`str`)

    :rtype: dict

    * Code originally ported from IRISMustangMetrics R Cran Package
    (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, March 20). CRAN-Package IRISMustangMetrics.
    The Comprehensive R Archive Network. Retrieved from
    https://cran.r-project.org/web/packages/IRISMustangMetrics/index.html) and augmented and adapted for use within
    Pycheron

    "The correlation is a value in the range [+/-1]. The 'pearson r' correlation is a measure of the strength
    and direction of the linear relationship between two variables that is defined as the (sample) covariance of the
    variables divided by the product of their (sample) standard deviations"
    (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf).
    Pearson's correlation requires that each trace (in this case) is normally distributed and continuous
    during the measurement. Positive correlations (+1) imply that as tr1 increases, so does tr2 (exact linear
    relationship). Negative correlations (-1) imply that as tr1 increases, tr2 decreases (exact linear relationship).
    While a value of 0 implies no correlation.
    (Last 4 sentences paraphrased from from scipy.stats.mstats.pearsonr:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mstats.pearsonr.html)

    "The p-value roughly indicates the probability of an uncorrelated system producing datasets that have a Pearson
    correlation at least as extreme as the one computed from these datasets"
    (https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mstats.pearsonr.html). The null hypothesis is that
    tr1 and tr2 are uncorrelated. The p-value is a number between zero and one and describes the statistical
    statistical significance of your observation. It represents the probability that your data would have arisen by
    random chance (i.e., null hypothesis is true) (http://www.eecs.qmul.ac.uk/~norman/blog_articles/p_values.pdf).
    Smaller p-values indicate that the null hypothesis should be rejected that your data would have arisen if the null
    hypothesis were true. Small values, closer to zero imply that the null hypothesis should be rejected. While larger
    values, imply that the null hypothesis is likly to be true. Common statistically significant p-values used are
    P < 0.05, P < 0.01, or P < 0.001.

    Seismic traces passed to `correlationMetric` should have the same network and station, span the same time range
    and have the same sampling rate (https://cran.r-project.org/web/packages/IRISMustangMetrics/IRISMustangMetrics.pdf).

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
        >>> {'p-value': 0.0, 'metric_name': 'correclationMetric', 'start_time': UTCDateTime(2013, 3, 1, 0, 0, 0, 19500),
             'correlation_coefficient': -0.2041723553359387, 'snclq': u'IU.ANMO.00.BH1.M:IU.ANMO.00.BH2.M',
             'end_time': UTCDateTime(2013, 3, 1, 23, 59, 59, 969500)}

    """

    # Set up logger
    if logger is None:
        logger = Logger(None)

    # -----------Compatibility Checks------------------------------------------
    # Start time and end time for trace data 1
    starttime = tr1.stats.starttime
    endtime = tr1.stats.endtime

    # Test if start times are the same for both traces
    if tr2.stats.starttime != starttime:
        logger.error("correlationMetric(): Start times do not match.")
        return

    # Test if end times are the same for both traces
    if tr2.stats.endtime != endtime:
        logger.error("correlationMetric(): End times do not match.")
        return

    # Test if sampling rates are the same for both traces
    if tr1.stats.sampling_rate != tr2.stats.sampling_rate:
        logger.error("correlationMetric(): Sampling Rates do not match.")
        return

    # Ensure traces are the same length
    l1 = len(tr1)  # Length of trace data 1
    l2 = len(tr2)  # Length of trace data 2

    # Complain if len of traces differs by more than 2 samples
    if abs(l2 - l1) > 2:
        logger.error("correlationMetric(): Incompatible lengths tr1 = %s, tr2 = %s" % (l1, l2))
        return
    # Get min length
    else:
        min_length = min(l1, l2)

    # Test if network and stations are the same for both traces
    if tr1.stats.network != tr2.stats.network or tr1.stats.station != tr2.stats.station:
        logger.error("correlationMetric(): Incompatible trace ids, tr1 = %s, tr2 = %s" % (tr1.get_id(), tr2.get_id()))
        return

    # Check for flat-lined data
    # Make a copy of tr1, tr2 as isDC applies tapering and demeaning techniques, etc that permanently change the trace
    # and we don't want to affect the correlation results
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

    # Merge snclq trace information. Station, network will be the same, as it's likely the mseed quality will be,
    # but decided to merge entire snclq information instead of just merging location and channel information,
    #  like R code does.
    # This has the same effect as the R code but just easier to take the full SNCL info with ObsPy
    try:
        snclq = (
            tr1.get_id() + "." + tr1.stats.mseed.dataquality + ":" + tr2.get_id() + "." + tr2.stats.mseed.dataquality
        )
    # If data quality is missing, just merge the sncl info
    except AttributeError:
        snclq = tr1.get_id() + ":" + tr2.get_id()

    # -----------Correlation Function------------------------------------------
    # Calculate the correlation metric and the P-value
    # Prior to the calculation ignore NaNs (if present), then feed into scipy's pearsonr function
    tr1 = tr1.data[~np.isnan(tr1.data)]
    tr2 = tr2.data[~np.isnan(tr2.data)]
    corr = pearsonr(tr1[:min_length], tr2[:min_length])
    # Output metric information
    c = {
        "snclq": snclq,
        "start_time": starttime,
        "end_time": endtime,
        "correlation_coefficient": corr[0],
        "p_value": corr[1],
        "metric_name": "correlationMetric",
    }

    # If database defined, insert metric information
    if database_config is not None:
        database = Database(**database_config)
        database.insert_metric(c)

    return c
