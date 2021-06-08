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

# Code base originally written by Steven Magana-Zook from Lawrence Livermore National Laboratories
# Augmented and adapted for use within Pycheron by Pycheron team

__all__ = ["FeatureGenerator"]

import numpy as np
from scipy.signal import welch
import math


class FeatureGenerator:
    def generate(self, trace):
        """

         Function to generate drop out fraction, drop out fraction, distinct values ratio, packet time
         bandwidth product, frequency sigma, and discontinuity max value.

         returns: returns tuple of generated features:

                 * dropout fraction is the number of discrete intervals in the trace with N or more consecutive samples
                   having the same value
                 * distinct values is the number of distinct values divided by the total number of values (this is
                   inversely related to the quantization error)
                 * tbp is packet time bandwidth product
                 * frequency sigma - standard deviation of the frequency
                 * discontinuity max value - maximum value of discontinuities

         rtype: tuple

         * Code base originally written by Steven Magana-Zook from Lawrence Livermore National Laboratories
           Augmented and adapted for use within Pycheron by Pycheron team

        ** Resources for ML/features selected Method **:

        #. Dodge, D.A., and D.B. Harris (2016). Local-Scale Test of Dynamic Correlation Processors: Implications for
        Correlation-Based Seismic Pipelines, Bull Seis. Am. 106(2), pp. 435-452.

        #. Dodge, D.A., and W.R. Walter (2015). Initial Global Seismic Cross-Correlation Results: Implications for
        Empirical Signal Detectors, Bull. Seis. Am. 105(1), pp.240-256.

        """

        # Generate dropout fraction, distinct values, tbp, frequency sigma, and discontinuity max value features
        drop_out_fraction = self.__get_drop_out_fraction(trace)
        distinct_vals = self.__get_distinct_values_ratio(trace)
        tbp, freq_sigma = self.__get_tbp(trace)
        disc_max_value = self.__get_disc_max_value(trace)

        # If tbp is nan, set it to 0
        if tbp == np.nan or math.isnan(tbp):
            tbp = 0
        # If frequency sigma is nan, set it to 0
        if freq_sigma == np.nan or math.isnan(freq_sigma):
            freq_sigma = 0

        # The order of return values MUST match the expected feature order of the model or the predictions will be
        # wrong. Return tuple of drop out fraction, distinct values, tbp, frequency sigma, and discontinuity max value.
        return drop_out_fraction, distinct_vals, tbp, freq_sigma, disc_max_value

    def __get_drop_out_fraction(self, trace, min_dropout_window_length=10):
        """

        Function to find the drop out fraction within a given trace. The dropout fraction is the number of discrete
        intervals in the trace with N or more consecutive samples having the same value

        :param tr: obspy trace object
        :type tr: obspy.core trace object
        :param min_dropout_window_length: threshold of min # of samples to tolerate in a sequence before considered
                                          a dropout. Thus discrete intervals are not included in the drop outs count
                                          if they are less than this minimal threshold
        :type min_dropout_window_length: int

        return: For each extreme value, the standard deviation from the rest of the trace is calculated. This
                returns the max standard deviation of the extreme values found in the trace.

        rtype: float

        """

        # Subtract each i+1 element from each ith element
        diffs = np.zeros(len(trace.data) - 1)
        differences = diffs[0 : len(diffs)] = (
            trace.data[1:] - trace.data[0 : len(trace.data) - 1]
        )

        # Initialize dropout counter at 0 and create sequence_lengths list to append to
        # Sequence length drop out counter = 1 when differences are 0, meaning they have the same value
        # Sequence length drop out counter = 0 when differences aren't 0, meaning they do not have the same value
        dropout_counter = 0
        sequence_lengths = []
        # Loop through differences and calculate
        for i in range(0, len(differences) - 1):
            if differences[i] == 0:
                # This i-th index was the same value as the previous entry
                dropout_counter += 1
            else:
                # New value - save off the previous series and reset counting
                sequence_lengths.append(dropout_counter)
                dropout_counter = 0

            # Edge case: At the end of the sequence
            if i == len(differences) - 2:
                sequence_lengths.append(dropout_counter)

        # Do not include in the drop outs counts if its less than our minimal threshold
        sequence_lengths = np.where(
            np.array(sequence_lengths) < min_dropout_window_length, 0, sequence_lengths
        )

        # Drop out fraction is the number of total counts that we considered drop outs as a ratio of the length of the
        # stream
        return np.sum(sequence_lengths) / trace.data.shape[0]

    def __get_distinct_values_ratio(self, trace):
        """
        Function to calculate the distinct values ratio of input trace. Distinct values is the number of distinct values
        divided by the total number of values (this is inversely related to the quantization error)

        :param tr: obspy trace object
        :type tr: obspy.core trace object

        returns: Number of distinct values divided by the total number of values
        rtype: float

        """
        # Calculate distinct values ratio, which is the number of distinct values divided by the total number of
        # values in the trace
        return np.unique(trace.data).shape[0] / trace.data.shape[0]

    def __get_tbp(self, trace):
        """
        Function to calculate the time bandwidth product (tbp) and frequency sigma of input trace.

        :param tr: obspy trace object
        :type tr: obspy.core trace object

        returns: tbp value and frequency sigma
        rtype: tuple

        """
        # Calculate the psd of the trace data, get the sampling frequency interval spacing
        f, Pxx = welch(trace.data, trace.stats.sampling_rate)
        sample_interval_freq = f[1] - f[0]

        # Calculate the spectral centroid and the frequency sigma using the sampling frequency interval spacing
        omega = self.__compute_centroid(Pxx, sample_interval_freq)
        sigmaf = self.__compute_sigma(Pxx, sample_interval_freq, omega)

        # Calculate the spectral centroid using the data sampling interval spacing.
        # Calculate sigma of the data
        t0 = self.__compute_centroid(Pxx, trace.stats.delta)
        sigma = self.__compute_sigma(trace.data, trace.stats.delta, t0)

        # Return tbp and frequency sigma
        return sigma * sigmaf, sigmaf

    def __get_disc_max_value(self, trace, factor=10, max_range=1e5):
        """
        Function to find the max discontinuity within a given trace object

        :param tr: obspy trace object
        :type tr: obspy.core trace object
        :param factor: scale factor used to scale the calculated differences standard deviation (default = 10) by to
                       create a threshold for determining extreme values. This value is used as the discontinuity
                       threshold.
        :type factor: int
        :param max_range: Maximum range for the trace's data
        :type max_range: float

        return: For each extreme value, the standard deviation from the rest of the trace is calculated. This
                returns the max standard deviation of the extreme values found in the trace.

        """

        # Calculate range of data
        mrange = np.max(trace.data) - np.min(trace.data)
        # If mrange is greater than max range, scale the trace data by max_range/mrange
        if mrange > max_range:
            scale = max_range / mrange
            trace.data = trace.data * scale

        # Subtract each i+1 element from each ith element then square those differences
        diffs = np.zeros(len(trace.data) - 1)
        differences = diffs[0 : len(diffs)] = (
            trace.data[1:] - trace.data[0 : len(trace.data) - 1]
        )
        differences = differences ** 2

        # Compute statistics on the differences to generate a discontinuity threshold
        diff_avg = np.mean(differences)
        diff_std = np.std(differences)
        threshold = factor * diff_std

        # Call any value an "extreme value" if it's value minus the mean is above the threshold
        extreme_values = differences[differences - diff_avg > threshold]
        if len(extreme_values) == 0:
            return 0

        # Find out how many standard deviations away each extreme value is from the rest of the trace.
        std_away = np.array(
            [0 if diff_std == 0 else num / diff_std for num in extreme_values]
        )

        # Return the max standard deviation
        return np.max(std_away)

    def __get_definite_integral(self, data, delta, idx1, idx2):
        """
        Calculate the definite integral of the input data
        :param data: Power spectral density of input trace data or input trace data (i.e., tr.stats.data)
        :type data: numpy.ndarray
        :param delta: sampling interval of frequency or data (i.e., f[1]-f[0] from Pxx frequency array or tr.stats.delta
        :type delta: float
        :param idx1: Lower integral limit
        :type idx1: int
        :param idx2: Upper integral limit
        :type idx2: int

        returns: returns the definite integral of the input data
        rtype: float

        """

        # Initialize result as data of lower integral limit index + data of upper integral limit index divided by two
        result = (data[idx1] + data[idx2]) / 2
        # Initialize j to be idx1 (integral lower limit) + 1
        j = idx1 + 1
        # While j is less than integral upper limit
        while j < idx2:
            # Add data[j] to the result
            result += data[j]
            # Increment j
            j += 1
        # Return final result * delta
        result *= delta
        return result

    def __compute_centroid(self, data, sample_interval):
        """
        Compute the spectral centroid using power spectral density and sampling interval. Sampling interval could
        be frequency sampling interval or data sampling interval (i.e., trace.stats.delta)

        :param data: Power spectral density of input trace
        :type data: numpy.ndarray
        :param sample_interval: sampling interval of frequency or data (i.e., f[1]-f[0] from Pxx frequency array, or
                                tr.stats.delta)
        :type sample_interval: float

        return: spectral centroid (raw centroid / energy)
        rtype: float

        """

        # Compute energy of input data based on sampling interval, using 0 as lower limit of integral and len(data) - 1
        # as upper limit of integral, if energy is 0, return nan
        energy = self.__get_definite_integral(
            np.abs(data), sample_interval, 0, len(data) - 1
        )
        if energy == 0:
            return np.nan

        # Calculate centroid kernel and raw centroid
        centroidKernel = (np.array(range(len(data))) * sample_interval) * np.abs(data)
        rawCentroid = self.__get_definite_integral(
            centroidKernel, sample_interval, 0, len(centroidKernel) - 1
        )

        # Centroid is the raw centroid over energy
        return rawCentroid / energy

    def __compute_sigma(self, data, delta, centroid):
        """
        Compute standard deviation (sigma) of input data

        :param data: Power spectral density of input trace data or input trace data (i.e., tr.stats.data)
        :type data: numpy.ndarray
        :param delta: sampling interval of frequency or data (i.e., f[1]-f[0] from Pxx frequency array or tr.stats.delta
        :type delta: float
        :param centroid: spectral centroid
        :type centroid: float

        return: standard deviation of input data
        rtype: float

        """
        # Treat data as absolute value of data through calculations
        data = np.abs(data)

        # Compute energy of input data based on sampling interval, using 0 as lower limit of integral and len(data) - 1
        # as upper limit of integral
        # If energy is 0, return nan
        energy = self.__get_definite_integral(data, delta, 0, len(data) - 1)
        if energy <= 0:
            return np.nan

        # Get variance kernel
        idxs = np.array(range(len(data)))
        tmt0 = (idxs * delta) - centroid
        varianceKernel = tmt0 * tmt0 * data

        # Compute variance
        variance = self.__get_definite_integral(
            varianceKernel, delta, 0, len(varianceKernel) - 1
        )
        variance = variance / energy

        # Return standard deviation
        return np.sqrt(variance)
