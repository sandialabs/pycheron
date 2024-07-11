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

__all__ = ["noiseModelInfra"]

import os
import numpy as np
from scipy.interpolate import interp1d


def noiseModelInfra(freq):
    """
    The noiseModelInfra function returns the IDC High Noise Model and IDC Low Noise Model from the 2012 Brown et. al. paper.
    Returned values are based on user specified frequency argument. [#]

    :param freq: array of frequencies
    :type freq: numpy.ndarray

    :return:

        * idclnm - IDC low noise model
        * idchnm - IDC high noise model

    :rtype:

        * numpy.ndarray
        * numpy.ndarray

    * IDC infrasound low and high noise models obtained from Brown et. al., 2012: 
    "Brown, D., Ceranna, L., Prior, M. et al. 
    The IDC Seismic, Hydroacoustic and Infrasound Global Low and High Noise Models. Pure Appl. Geophys. 171, 361–375 (2014)". 
    
    **Example**

    .. code-block:: python

        from pycheron.psd.noise.noiseModelInfra import noiseModelInfra
        from pycheron.psd.noise.getNoise import getNoise
        from pycheron.psd.psdList import psdList

        # TODO update test data to include an infrasound file; example IM.I53H1.BDF to compare to
        #test data
        data = 'test/test_data/6e_sp06_ehe.407438.tar.mseed'

        # reading in stream
        st = obspy.read(data)

        # calculating psds
        psds = psdList(st)
        # Get instrument corrected psds
        f,n,psd = getNoise(psds, units='DEF')

        freq = f[0]
        period = 1/freq
        # calculating idchnm, idclnm
        idclnm,idchnm = noiseModelInfra(freq)

    **Plotting**

    .. code-block:: python

        import matplotlib.pyplot as plt

        # Normally shown on log xscale need to update
        plt.plot(idclnm[::-1])
        plt.plot(idchnm[::-1])
        plt.title('IDC Noise Models')
        plt.xlabel('Period (s)')
        plt.ylabel('Pa^2/Hz (dB)')


    .. image:: _static/noiseModelIDC.png

    """

    # Get directory path for loading IDC low/high noise models
    dir_path = os.path.dirname(os.path.realpath(__file__))
 
    # Load files with noise models & frequency. Column 1 is log10(Freq), Column 2 is log10 PSD (low noise/high noise model Pa^2/Hz)
    lnm = np.loadtxt(dir_path + '/noise_models/IDC2012_Low.txt', skiprows = 1)
    hnm = np.loadtxt(dir_path + '/noise_models/IDC2012_High.txt', skiprows = 1)
    
    # Create interpolated function for low and high noise models for use in interpolating and extrapolating (if needed) frequency input. Remove log10 from frequency values and mulitply noise models by 10 to get into right dB units. 
    low_noise_interp = interp1d(10**lnm[:,0], 10*lnm[:,1], fill_value = 'extrapolate')
    high_noise_interp = interp1d(10**hnm[:,0], 10*hnm[:,1], fill_value='extrapolate')

    # Create interpolated low/high noise models based on user input frequency and return the values
    idclnm = low_noise_interp(freq)
    idchnm = high_noise_interp(freq)

    return idclnm, idchnm
