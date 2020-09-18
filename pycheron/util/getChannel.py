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

__all__ = ["getChannel", "filterChannel", "getChannelName"]


def getChannel(st, channel):
    """

    Gets channels from stream with multiple channels

    :param st: Obspy stream object
    :type st: obspy.core.stream.Stream
    :param channel: Name of channel component - E, N, Z, or component channel name - EHE, EHN, EHZ, LHE, LHN,
                    LHZ, LOG, VM1, VM2, VM3, etc. if you want a specific component. Alternatively can specify
                    EH\* if want to select all components with a common band/instrument code.
    :type channel: str

    :return: Returns a new stream object containing the specific trace(s) with the specified channel.
    :rtype: obspy.core.stream.Stream

    **Example**

    .. code-block:: python

        from pycheron.util.getChannel import getChannel

        # Grab test data
        data = 'test/test_data/6e_sp06_Eall.397679.tar.mseed'
        data1 = 'test/test_data/7A_CABN_ALL.988887.tar.mseed'

        # Reading in stream then merge together (only merging for better example stream (e.g., stream with multiple
        # sncls))
        st = obspy.read(data)
        st1 = obspy.read(data1)
        ST = st+st1

        # getting EHE channel trace
        tr = getChannel(ST, 'EHE')
        # Output
        # 2 Trace(s) in Stream:
        # 6E.SP06..EHE | 2014-01-19T00:00:00.005000Z - 2014-01-19T21:01:07.985000Z | 100.0 Hz, 7566799 samples
        # 6E.SP06..EHE | 2014-01-19T21:07:56.000000Z - 2014-01-20T00:00:00.990000Z | 100.0 Hz, 1032500 samples

        # getting all E component
        tr = getChannel(ST,'E')
        # Output
        # 3 Trace(s) in Stream:
        # 6E.SP06..EHE | 2014-01-19T00:00:00.005000Z - 2014-01-19T21:01:07.985000Z | 100.0 Hz, 7566799 samples
        # 6E.SP06..EHE | 2014-01-19T21:07:56.000000Z - 2014-01-20T00:00:00.990000Z | 100.0 Hz, 1032500 samples
        # 7A.CABN..BHE | 2013-11-01T00:00:00.000000Z - 2013-11-01T23:59:59.975000Z | 40.0 Hz, 3456000 samples

    """
    # Grab out specified component within stream object

    # First option is for specifying want all E component data, this could simply be EHE, or if multiple E component
    # channels could mean grab all E component data (e.g., EHE, LHE, etc)
    if len(channel) == 1:
        traces = st.select(component=channel)
        return traces
    # If option specified is a specific channel code,e.g., EHE or LHE, then only grab out that one. Alternatively, grab
    # out all components with a common band/instrument code (e.g., EH*, etc.)
    else:
        traces = st.select(id=str("*") + channel)
        return traces


def getChannelName(snclq):
    """
    Extracts the channel name from a SNCLQ id

    :param snclq: snclq id
    :type snclq: str

    :return: The name of the channel
    :rtype: str
    """
    chanName = snclq.split(".")[3]
    return chanName


def filterChannel(st):
    """
    Filters SOH channel traces from stream (ACE, LCE, LCQ, LOG, OCF, VCO, VEA, VEC, VEP, VKI, VMU, VM1, VMZ, VMV, VM2,
    VMN, VMW, VM3, VME, VPB)

    :param st: Obspy stream object
    :type st: obspy.core.stream.Stream

    :return: New stream with SOH channels removed.
    :rtype: obspy.core.stream.Stream
    """
    soh_chan = [
        "ACE",
        "LCE",
        "LCQ",
        "LOG",
        "OCF",
        "VCO",
        "VEA",
        "VEC",
        "VEP",
        "VKI",
        "VMU",
        "VM1",
        "VMZ",
        "VMV",
        "VM2",
        "VMN",
        "VMW",
        "VM3",
        "VME",
        "VPB",
    ]
    stN = st.copy()
    for i in range(len(st)):
        tr = st[i]
        snclq = tr.get_id()
        chan = getChannelName(snclq)
        if chan in soh_chan:
            stN.remove(tr)
    return stN
