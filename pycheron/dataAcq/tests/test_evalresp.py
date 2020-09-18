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

from obspy.core.utcdatetime import UTCDateTime
from ispaq.evalresp import evalresp as read_Evalresp


def test_read_Evalresp():
    vals = {
        "network": "UU",
        "station": "JLU",
        "location": "",
        "channel": "HHZ",
        "starttime": UTCDateTime("2006,304,00:00:00"),
        "minFreq": 0.005255602595335718,
        "maxFreq": 46.95060701207919,
        "nFreq": 106,
        "units": "VEL",
        "filename": "test_data/RESP.UU.JLU..HHZ",
    }

    iris = read_Evalresp(
        sfft=vals["minFreq"],
        efft=vals["maxFreq"],
        nfft=vals["nFreq"],
        filename=vals["filename"],
        date=vals["starttime"],
        station=vals["station"],
        channel=vals["channel"],
        network=vals["network"],
        locid=vals["location"],
        units=vals["units"],
        debug=True,
    )

    assert isinstance(iris, tuple)
    assert len(iris) > 0
