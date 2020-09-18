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

import numpy as np
from scipy import signal
from numpy.lib import stride_tricks


def sta_lta(wav, fs, t_short=1, t_long=10, method="squared"):
    """
    Short Time Average / Long Time Average (STA/LTA)
    :param wav: (array like)
    :param fs: sampling frequency
    :param t_short: duration of short time window (seconds)
    :param t_long: duration of long time window (seconds)
    :param method: either name (str) or index (int)
    :return: (array like) cropped by a total of (t_short + t_long)
                reason being that an initial (t_short + t_long) seconds is
                needed to compute the first STA/LTA value
    """
    methods = [
        "absolute",
        "squared",
        "difference",
        "modified energy ratio",  # very sharp
        "allen",  # smooth 'squared'
        "BCTF",
    ]

    if isinstance(method, int):
        assert method < len(
            methods
        ), "not a valid method for sta_lta(), " "max = " + str(len(methods) - 1)
        method = methods[method]
    assert method in methods, "'" + method + "' is not a valid method " "for sta_lta()"

    # pad signal with random samples to prevent large spike at beginning
    # x = np.append(rms(wav) * np.random.randn(int(t_long * fs)), wav)
    # x = append(zeros(int(t_long * fs)), wav)
    x = wav

    if method == "absolute":
        x = abs(x)
    elif method == "allen":
        b = [1, -1]
        x_diff = signal.lfilter(b, 1, x)
        x = x ** 2 + 3 * x_diff ** 2
    elif method == "BCTF":
        b = [1, -1]
        x_diff = signal.lfilter(b, 1, x)

        envelope = x ** 2 + x_diff ** 2
        envelope *= np.cumsum(x ** 2) / np.cumsum(x_diff ** 2)
        x = envelope ** 2

        x_mean = np.zeros(len(x))
        x_var = np.zeros(len(x))

        for i in range(0, len(x)):
            x_mean[i] = np.mean(x[0 : (i + 1)])
            x_var[i] = np.var(x[0 : (i + 1)])

        x = (x - x_mean) / x_var
    else:
        x **= 2

    sta_win_len = int(t_short * fs)
    lta_win_len = int(t_long * fs)

    b = np.ones(sta_win_len) / sta_win_len
    sta = signal.lfilter(b, 1, x)

    b = np.ones(lta_win_len) / lta_win_len
    lta = signal.lfilter(b, 1, x)

    # shift sta ahead of lta
    sta = sta[sta_win_len::]
    lta = lta[0:-sta_win_len]

    sta_lta_val = sta / lta

    if method == "difference":
        b = [1, -1]
        sta_lta_val = signal.lfilter(b, 1, sta_lta_val)

    elif method == "modified energy ratio":
        sta_lta_val = (abs(wav) * sta_lta_val) ** 3

    return sta_lta_val[int(t_long * fs) : :]


def stft(sig, frame_size, fs=1, overlap_factor=0.5, window=np.hanning):
    """
    Short Time Fourier Transform (STFT) Function
    :param sig: input waveform
    :param frame_size: stft frame length (samples)
    :param fs: sampling frequency
    :param overlap_factor: (< 1)
    :param window: window type
    :return:
      X         - stft X[n,k]
      t_stft    - time vector
      f_stft    - frequency vector
    """
    frame_size = int(frame_size)

    win = window(frame_size)
    win /= sum(win)
    hop_size = int(frame_size - np.floor(overlap_factor * frame_size))

    # zeros at beginning (thus center of 1st window should be for sample nr. 0)
    samples = np.append(np.zeros(int(np.floor(frame_size / 2.0))), sig)
    # cols for windowing
    cols = int(np.ceil((len(samples) - frame_size) / float(hop_size)) + 1)
    # zeros at end (thus samples can be fully covered by frames)
    samples = np.append(samples, np.zeros(int(frame_size)))

    frames = stride_tricks.as_strided(
        samples,
        shape=(int(cols), int(frame_size)),
        strides=(samples.strides[0] * hop_size, samples.strides[0]),
    ).copy()
    frames *= win

    t_stft = np.arange(0, frames.shape[0]) * hop_size / fs
    # fft.rfft returns only positive frequencies
    f_stft = np.arange(0, np.floor(frames.shape[1] / 2) + 1) / frames.shape[1] * fs

    X = np.fft.rfft(frames)

    return X, t_stft, f_stft
