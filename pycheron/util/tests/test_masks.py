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

import pytest
import numpy as np
from obspy.core.utcdatetime import UTCDateTime
from pycheron.util.masks import consecutive, samples2time


class MockTraceVals:
    def __init__(self, fs, index, start):
        self.fs = fs
        self.index = index
        self.start = start


def test_consecutive():
    # TODO: Consider refactoring the function this tests to something simpler
    # All it should be doing is breaking each individual element in a np.array to be
    # Check deadChanMeanMetric.py for example
    # its own np array (ex: [1,5,11,17] -> [array([1]), array([5]), array([11]), array([17])] )
    # There should be a simpler way to implement this
    dat_array = np.array([1, 5, 11, 17])

    consec_vals = consecutive(dat_array)

    assert dat_array.size == len(consec_vals)

    for np_arr in consec_vals:
        assert np_arr.size == 1
        assert isinstance(np_arr, np.ndarray)


def test_samples2time(file_assets):
    mockTraceVals = MockTraceVals(
        fs=100.0,
        index=np.arange(1, 1000),
        start=UTCDateTime(2011, 1, 1, 2, 16, 8, 968300),
    )
    mask = samples2time(mockTraceVals.index, mockTraceVals.fs, mockTraceVals.start)

    assert len(mask) == len(mockTraceVals.index)
    assert isinstance(mask[0], dict)
    assert list(mask[0].keys()) == ["start_time", "end_time"]

    try:
        UTCDateTime(mask[0]["start_time"])
        UTCDateTime(mask[0]["end_time"])
    except Exception as e:
        pytest.fail("Dict values not returning datetime string: {e}".format(e=e))
