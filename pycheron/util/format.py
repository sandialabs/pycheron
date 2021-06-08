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


def parse_snclq(snclq, comb=False):
    """
    This function takes a sncql string and parses it into individual variables.

    :param sncql: Formatted ex. IU.ANMO.00.BHZ.M
    :type snclq: str
    :param comb: If inputing 2 snclqs from correlation metrics formatted ex. SNCLQ1:SNCLQ2
    :type comb: bool

    :return: network, station, channel, location, quality. If ``comb=True``, will return 2 formatted SNCLQs to be run
             again individually.
    :rtype: str, str, str, str, str

    **Example**

    1. Individual SNCLQ

    .. code-block:: python

        from pycheron.util.format import parse_snclq
        snclq = "7A.CABN..BHE"
        network, station, channel, location, quality = parse_snclq(snclq)

    2. Combined SNCLQ

    .. code-block:: python

        from pycheron.util.format import parse_snclq
        snclq_comb = "7A.CABN..BHE:7A.CABN..BHN"
        snclq1, snclq2 = parse_snclq(snclq_comb)
        network1, station1, channel1, location1, quality1 = parse_snclq(snclq1)
        network2, station2, channel2, location2, quality2 = parse_snclq(snclq2)

    """
    # If comb, parse the sncl into two sncl objects and return them
    # If inputing 2 snclqs from correlation metrics formatted ex. SNCLQ1:SNCLQ2
    if comb:
        parse = snclq.split(":")
        snclq1 = parse[0]
        snclq2 = parse[1]

        return snclq1, snclq2
    # Otherwise just parse out the network, station, channel, location, and quality information
    else:
        parse = snclq.split(".")
        network = parse[0]
        station = parse[1]
        location = parse[2]
        channel = parse[3]
        quality = ""
        if len(parse) == 5:
            quality = parse[4]

        return network, station, channel, location, quality


def format_mseed_filenames(data_dir, new_data_dir):
    """
    This function takes in a data directory path with mseed files and reads in each file to reformat each stream into
    the correct format ``callPycheronMetric`` needs. (i.e. each filename is ``<station>_<channel>_<julian day>.mseed``).
    Then outputs them to a new data folder.

    :param data_dir: Path to the current data directory contain all the mseed files
    :type data_dir: str
    :param new_data_dir: Name of directory to save newly formatted files.
    :type new_data_dir: str

    **Example**

    1. Combined Channels and Days

       This example shows how to use this function on a folder containing mseed files with combined channels/dates.

       *Current Folder Structure*

       .. code-block:: console

           # data folder containing a singular mseed file or a singular station
           data_dir
               ANMO_All_chans_3_days.mseed

       *New Folder Structure*

       .. code-block:: console

           new_data_dir
                ANMO_EHE_001.mseed
                ANMO_EHN_001.mseed
                ANMO_EHZ_001.mseed
                ANMO_EHE_002.mseed
                ANMO_EHN_002.mseed
                ANMO_EHZ_002.mseed
                ANMO_EHE_003.mseed
                ANMO_EHN_003.mseed
                ANMO_EHZ_003.mseed

    2. Nested Folder

       This example shows how to use this function on a nested folder.

       *Current Folder Structure*

       .. code-block:: console

           # datafolder containing station folders
           data_dir
                # Station folder containing individual channel folders with individual day .mseed files
                BGU
                    EHE
                        day1.mseed
                        day2.mseed
                    EHN
                        day1.mseed
                        day2.mseed
                    EHZ
                        day1.mseed
                        day2.mseed
                # Station folder containing individual channels, but combined days
                JPU
                    BHN_all_days.mseed
                    BHE_all_days.mseed
                    BHZ_all_days.mseed
                # Station folders containing individual days, but combined channels
                TCRU
                    EH*_day1.mseed
                    EH*_day2.mseed

       *New Folder Structure*

       .. code-block:: console

           new_data_dir
                BGU_EHE_001.mseed
                BGU_EHN_001.mseed
                BGU_EHZ_001.mseed
                BGU_EHE_002.mseed
                BGU_EHN_002.mseed
                BGU_EHZ_002.mseed
                JPU_EHE_001.mseed
                JPU_EHN_001.mseed
                JPU_EHZ_001.mseed
                JPU_EHE_002.mseed
                JPU_EHN_002.mseed
                JPU_EHZ_002.mseed
                TCRU_EHE_001.mseed
                TCRU_EHN_001.mseed
                TCRU_EHZ_001.mseed
                TCRU_EHE_002.mseed
                TCRU_EHN_002.mseed
                TCRU_EHZ_002.mseed

    """
