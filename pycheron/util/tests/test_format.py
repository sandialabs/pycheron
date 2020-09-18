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
from pycheron.util.format import parse_snclq


def test_parse_snclq():
    # Test combination set to False
    snclq = "7A.CABN..BHE"
    network, station, channel, location, quality = parse_snclq(snclq, comb=False)

    assert network == "7A"
    assert station == "CABN"
    assert channel == "BHE"
    assert location == ""
    assert quality == ""

    # Test combination set to True
    snclq_comb = "7A.CABN..BHE:7A.CABN..BHN"
    snclq1, snclq2 = parse_snclq(snclq_comb, comb=True)
    network1, station1, channel1, location1, quality1 = parse_snclq(snclq1)
    network2, station2, channel2, location2, quality2 = parse_snclq(snclq2)

    assert network1 == "7A"
    assert station1 == "CABN"
    assert channel1 == "BHE"
    assert location1 == ""
    assert quality1 == ""

    assert network2 == "7A"
    assert station2 == "CABN"
    assert channel2 == "BHN"
    assert location2 == ""
    assert quality2 == ""
