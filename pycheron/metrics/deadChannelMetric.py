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

from pycheron.util.logger import Logger
from pycheron.psd.noise.deadChannel import isDC

__all__ = ["deadChannelMetric"]


def deadChannelMetric(st, logger=None, database=None):
    """
    Individaul metric that tests a channel to see if it is a dead channel

    :param st: Stream object
    :type st: obspy.core.stream.Stream
    :param logger: logger object
    :type logger: pycheron.util.logger.Logger
    :param database: database object
    :type database: pycheron.db.sqllite_db.Database

    :return: list of dictionaries for each channel in stream with the following keys and types

                * is_dead_channel (`bool`)
                * snclq (`str`)
                * start_time (`str`)
                * end_time (`str`)
                * metric_name (`str`)

    :rtype: list

    **Examples**

    .. code-block:: python

        # Initialize IRIS client
        client = Client("IRIS")

        # Define start/end time and then get stream object
        starttime = UTCDateTime("2012-12-12T00:00:00.000")
        endtime = UTCDateTime("2012-12-13T00:00:00.000")
        st = client.get_waveforms("AK","GHO","","BH*",starttime, endtime)

        dc = deadChannelMetric(st)
        for i in range(len(dc)):
            print "is_dead_channel: ", dc[i]["is_dead_channel"]
            print "sncql: ", dc[i]['snclq']
            print "start_time: " , dc[i]['start_time']
            print "end_time: ", dc[i]['end_time']
            print "metric_name: ", dc[i]['metric_name']
            print "---------"

        >>> is_dead_channel:  False
        >>> sncql:  AK.GHO..BHE
        >>> start_time:  2012-12-12T00:00:00.008400
        >>> end_time:  2012-12-12T23:59:59.988400
        >>> metric_name:  deadChannelMetric
        >>> ---------
        >>> is_dead_channel:  False
        >>> sncql:  AK.GHO..BHN
        >>> start_time:  2012-12-12T00:00:00.008400
        >>> end_time:  2012-12-12T23:59:59.988400
        >>> metric_name:  deadChannelMetric
        >>> ---------
        >>> is_dead_channel:  False
        >>> sncql:  AK.GHO..BHZ
        >>> start_time:  2012-12-12T00:00:00.008400
        >>> end_time:  2012-12-12T23:59:59.988400
        >>> metric_name:  deadChannelMetric
        >>> ---------

    """
    if logger == None:
        logger = Logger(None)

    all = []
    for i in range(len(st)):
        tr = st[i]
        dc = isDC(tr, False, True, 0)
        data = {
            "is_dead_channel": dc,
            "snclq": tr.get_id(),
            "start_time": tr.stats.starttime.isoformat(),
            "end_time": tr.stats.endtime.isoformat(),
            "metric_name": "deadChannelMetric",
        }
        all.append(data)

    if database is not None:
        database.insert_metric(all)
    return all
