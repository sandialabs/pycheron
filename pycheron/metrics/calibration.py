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

# TODO this function needs a major rewrite, it's too much inside one function and could be better modularized
# TODO need to look into threshold parametrizations
# TODO need masking for cals


import os
import sys

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path + "/..")
import obspy
import numpy as np

from pycheron.db.sqllite_db import Database

# from scipy import signal
import pickle
import warnings
import joblib
from pycheron.metrics.utils import *
from pycheron.metrics.calibration_plotting_utils import *


warnings.filterwarnings("ignore")

__all__ = [
    "LogComressor",
    "leaky_integrator",
    "LTIFilter",
    "FrameContainer",
    "MLApplier",
    "CalibrationHelper",
    "unscreened_calibration",
    "calibration",
    "calibrationMetric",
]

__calibration_demo__ = False


class CustomUnpickler(pickle.Unpickler, object):
    def find_class(self, module, name):
        if module == "__main__":
            module = "pycheron.metrics.calibration"
        return super(CustomUnpickler, self).find_class(module, name)


class LogComressor(object):
    """Class to log compress data, then replace all instances of -INF with minimum of
    that row, and replace INF with maximum of that row. Class used to create compressor.pkl"""

    @staticmethod
    def _compress(features):
        with np.errstate(divide="ignore"):
            features = np.log(features.copy())
        return features

    def fit(self, features):
        features = self._compress(features)

        self._mins_and_maxs = []

        for j in range(0, features.shape[1]):
            temp = features[:, j]
            temp = temp[abs(temp) != np.inf]
            self._mins_and_maxs.append([np.min(temp), np.max(temp)])

        return self

    def transform(self, features):
        assert hasattr(self, "_mins_and_maxs"), (
            "LogCompressor object has " "not been fit yet!"
        )
        assert len(self._mins_and_maxs) == features.shape[1], "invalid input sizes"

        features = self._compress(features)

        for i in range(0, features.shape[0]):
            for j in range(0, features.shape[1]):
                if features[i, j] == np.inf:
                    features[i, j] = self._mins_and_maxs[j][1]
                elif features[i, j] == -np.inf:
                    features[i, j] = self._mins_and_maxs[j][0]

        return features


def leaky_integrator(fs, dur):
    """
    Single pole filter (leaky integrator I think)
    The impulse response is a decaying exponential.  The longer the <dur>
    is, the longer this exponential will take to decay.
    :param fs: sampling frequency (Hz)
    :param dur: time at which impulse response has decayed to below 90%
                longer <dur> implies longer startup and more aggressive LPF
    :return: filter coefficients
    """
    B = [20 / fs * 0.01150 * 10 / dur]
    A = [1, B[0] - 1]

    return B, A


class LTIFilter(object):
    """
    Helper class that manages storing and updating initial/final conditions.
    """

    def __init__(self, A, B, zi=None):
        self.A = A
        self.B = B
        self.zi = zi

    def has_initial_condition(self):
        return self.zi is not None

    def filter(self, x):
        y, zf = signal.lfilter(self.B, self.A, x, zi=self.zi)

        self.zi = zf

        return y


class FrameContainer(object):
    @property
    def max_size(self):
        return self._max_size

    def __init__(self, max_size):
        self._max_size = int(max_size)

        self.clear()

    def clear(self):
        self._data = []

    def size(self):
        return len(self._data)

    def empty(self):
        return self.size() == 0

    def full(self):
        return self.size() == self._max_size

    def add(self, item):
        if self.full():
            del self._data[0]

        self._data.append(item)

    def get(self, index):
        assert not self.empty(), "container is empty"
        assert self.size() > index, "index is higher than container size"

        return self._data[-1 - index]

    def __iter__(self):
        return iter(self._data[::-1])


class MLApplier(object):
    """This class screens detections using a trained machine learning
    classifier and some preprocessing steps that are all loaded from pickle
    files."""

    def __init__(self, loadpath):
        self.load_ml_objects(loadpath)

    def load_ml_objects(self, dirpath):
        filename = dirpath + "compressor.pkl"
        with open(filename, "rb") as fid:
            unpickler = CustomUnpickler(fid)
            compressor = unpickler.load()

        filename = dirpath + "imputer.joblib"
        imputer = joblib.load(filename)

        filename = dirpath + "scaler.joblib"
        scaler = joblib.load(filename)

        filename = dirpath + "classifier.joblib"
        classifier = joblib.load(filename)

        self.compressor = compressor
        self.imputer = imputer
        self.scaler = scaler
        self.classifier = classifier

    def iscal(self, detection):
        """Return True if it is a calibration"""

        # get data
        issue_start_time = float(detection["cal_start_time"])
        issue_end_time = float(detection["cal_end_time"])
        original_start_time = float(detection["initial_cal_start_time"])
        original_end_time = float(detection["initial_cal_end_time"])
        cal_duration = issue_end_time - issue_start_time
        Ts_minus_Tos = issue_start_time - original_start_time
        Toe_minus_Te = -issue_end_time + original_end_time

        features = [cal_duration, Ts_minus_Tos, Toe_minus_Te]
        features += [
            detection["max_ySmooth"],
            detection["max_stalta"],
            detection["env_mean"],
            detection["env_median"],
            detection["env_var"],
            detection["stft_var"],
        ]

        # transform features
        features = np.array(features).reshape(1, -1)
        features = self.compressor.transform(features)
        features = self.imputer.transform(features)
        features = self.scaler.transform(features)

        # apply classifier
        prediction = self.classifier.predict(features)

        return prediction == 1


class CalibrationHelper(object):
    """
    This class is a helper to the calibration function.  A single
    private object is instantiated and used by the function for
    configuration of parameters and storing data from previously viewed
    traces, in the case that data is being streamed to the function.

    The function can be dramatically changed by adjusting parameter values
    in _config().

    update_state() handles all stored data and filters based on whether the
    new trace is a contiguous block after the previously viewed trace.
    """

    def __init__(self):
        self.fs = None

        self.do_plot = False

        self.fc_snclq = None
        self.fc_starttime = -99
        self.fc_endtime = -99

        self._config()

    def update_state(self, trace):
        """Update anything that changed, such as sampling rate
        and check whether caches need to be cleared, which occurs when
        non-continuous traces are provided, i.e. different SNCLQ or
        discontiguous start/end times.
        """
        if self.fs != trace.stats.sampling_rate:
            self.fs = trace.stats.sampling_rate
            self._update_filters()

        self.bufferdur = trace.stats.endtime - trace.stats.starttime

        snclq = trace.get_id()
        if (
            self.fc_snclq != snclq
            or abs(self.fc_endtime - trace.stats.starttime) > self.fc_time_eps
        ):
            self.fc_snclq = snclq
            self.fc_starttime = -99
            self.fc_endtime = -99
            self._update_filters()
            self._update_frame_containers()

        try:
            assert (
                self.bufferdur > self.minDetectionDuration
            ), "each frame must last longer than the minimum detection duration"
        except AssertionError:
            return

    def _config(self):
        """
        Configure parameter settings for the calibration function.  The
        current configuration is tuned pretty well for IMS stations.  It is
        unknown if this performs well on other stations, however,
        the objective is to measure relative levels, i.e. how large is this
        value compared to typical values on this channel? Thus, it should
        adapt to any network.

        Any "dur..." for a filter is specifying the <dur> to be used in the
        design of a leaky_integrator filter.
        """
        # TODO (maybe) sweep for optimal ST, MT, LT and threshold
        # make it something that can be automatically done in the future

        # duration for high-pass filter of x
        # (high-pass is x - xLowPass)
        self.durXHP = 3600 * 24
        # self.durXHP = 60

        # short-time, long-time and threshold for LTA/STA ratio
        self.durShort = 15
        # self.durShort = 3
        self.durLong = 60
        self.threshold = 3  # TODO: Tested some other thresholds, like 2.5, to try to catch other calibrations
        # TODO: need to play with this more

        # filtering LTA/STA to smooth spikes
        self.durYSmooth = 15

        # window duration for envelope filter (hanning window)
        self.durEnvFilter = 10
        # minimum gap between detection start/end times
        self.minDetectionDuration = 20

        # this filter estimates typical envelope energy.
        # it is only updated per frame if no detections are found
        self.durAbsXLP = 60 * 60 * 6
        # lag for estimate of typical envelope energy such that if a
        # cal pulse occurs it does not unnaturally raise the denominator
        self.envelope_ratio_denominator_lag = 3600 * 24

        # STA and LTA times to recalibrate the start time
        self.recal_tLong = 30
        self.recal_tShort = 0.5

        # frame container durations (max amount of previous information stored)
        self.fc_dur = 3600 * 24
        # maximum gap between frames to not clear cache (seconds)
        self.fc_time_eps = 60

        assert self.minDetectionDuration > self.durEnvFilter, (
            "gap between start and end times for a detection must be "
            "greater than the envelope filter start-up"
        )

        # path to pickled machine learning objects used to screen initial
        # calibration detections
        self.ml_loadpath = os.path.dirname(__file__) + "/calibration_ml_objects/"

    def _update_filters(self):
        self.filters = {}

        # filter DC from signal
        B, A = leaky_integrator(self.fs, self.durXHP)
        self.filters["DC_remover"] = LTIFilter(A, B)

        # compute STA, LTA on |x[n]|
        B, A = leaky_integrator(self.fs, self.durShort)
        self.filters["sta"] = LTIFilter(A, B)

        B, A = leaky_integrator(self.fs, self.durLong)
        self.filters["lta"] = LTIFilter(A, B)

        # smooth LTA/STA
        B, A = leaky_integrator(self.fs, self.durYSmooth)
        self.filters["ltasta_smoother"] = LTIFilter(A, B)

        # computing signal envelope
        B = np.hanning(int(self.durEnvFilter * self.fs))
        B /= sum(B)
        self.filters["envelope"] = LTIFilter([1], B)

        # this filter estimates typical envelope level.
        # it is only updated per frame if no detections are found
        B, A = leaky_integrator(self.fs, self.durAbsXLP)
        self.filters["absxLP"] = LTIFilter(A, B)

    def _update_frame_containers(self):
        assert hasattr(self, "bufferdur"), "bufferdur hasn't been set yet!"

        # stores previous frames of x
        max_num_frames = self.fc_dur / self.bufferdur
        max_num_frames = max(int(max_num_frames), 1)

        self.FCx = FrameContainer(max_num_frames)
        self.FCy = FrameContainer(max_num_frames)

        self.FCsta = FrameContainer(max_num_frames)
        self.FClta = FrameContainer(max_num_frames)
        self.FCySmooth = FrameContainer(max_num_frames)

        max_num_frames = self.envelope_ratio_denominator_lag / self.bufferdur
        max_num_frames = max(int(max_num_frames), 1)
        self.FCzi = FrameContainer(max_num_frames)

    def iscal(self, detection):
        if not hasattr(self, "ml_applier"):
            self.ml_applier = MLApplier(self.ml_loadpath)

        return self.ml_applier.iscal(detection)


_chelper = CalibrationHelper()


def unscreened_calibration(trace, metric_store=None, debug=False, **params):
    """
    This function automatically detects calibrations in seismic data.  It was
    originally tuned to identify "stationary rect calibrations" which means
    high amplitude calibrations that have roughly a rect (square) envelope and
    appear stationary within the envelope. Upon running the code on data it
    was found that many other types of calibrations were also being identified.
    Thus, this algorithm is primarily tuned to a specific type of calibration
    but still has use in other cases.

    The general algorithm is as follows:
        1) an LTA/STA is run looking for triggers above threshold.  This is
        due to the fact that seismic arrivals can have abrupt transient
        ONSETs however the...OFFSET? is not as abrupt, whereas these
        calibrations to have an abrupt drop in energy.

        2) For these candidate triggers, start time is found by:
            a) computing an envelope for the data
            b) normalizing the envelope by the typical channel energy,
            such that anything above 1 is high
            c) searching BACKWARDS in time until something below 1 is
            identified, i.e. the channel has returned to normal energy levels

        3) The start and end times are refined.  The start time is refined
        by using an STA/LTA trigger, the end time is shifted to the time
        where LTA becomes > STA.  Both of these refined times will shorten
        the overall duration of the calibration

        4) These calibrations (could be 0 or more) can very often be false
        alarms. A second step is then implemented to screen the candidate
        calibrations.  The current implementation is somewhat complex,
        using a machine learning classifier and a collection of statistics,
        however, in the future it is possible to implement a rather simple
        set of conditions to keep a calibration.  For example:
            if duration < 5 seconds or energy < 3 X typical_energy_level:
                discard calibration

        5) The previous step is the final step to return a result, however
        some values are cached in order to process streams of data more
        accurately.  Details on this will be described shortly.

    Despite being a function, a private calibration helper object exists in
    volatile memory in order to assist calibration() as if it were an object.
    This object contains many parameter settings as well as caching data.
    One of the common reasons for caching data is that this algorithm
    identifies the END of a calibration, however, some calibrations can last
    very long.  If the calibration started before the provided trace's start
    time, more data is required to identify an accurate start time.

    ---------------------------------------------------------------------------
    About caching as using this function on sequential blocks of data.

    This function works best when applied to continuous streams of data.
    For example:

     snclq = <a unique and consistent snclq>
     for start_hour in range(0, 48, 1):
        end_hour = start_hour + 1
        trace = <a trace that starts at start_hour and ends at end_hour>
        results = calibration(trace)

     It is not so much important that the durations are constant, but that
     the end time of the previous frame is the same as the start tiem of the
     current frame.  If this is true, all cached data will be stored,
     not cleared, and can be used to assist in the calibration detection. This
     function was initially written to be processed on weeks to years at a
     time, however I would say at a maximum (so long as the current parameter
     settings are preserved) it takes a day to get all of the caches and data
     to steady-state.  Thus, if you want accurate results for Monday,
     run Sunday and Monday sequentially.

     The caches are FrameContainer objects, that much like a Queue that
     stores a fixed number of previous frames.  The important frame
     containers are:
        FCzi - this is the initial/final state of a filter used to estimate
        the typical (average) value of the signal envelope. This way the
        envelope can be normalized by this value and we can look for any
        times that the normalized envelope is significantly above one.

        FCx - stored time-series data (with DC removed). This is important
        for stepping back in time to find the start of a calibration,
        (in the event that the calibration spans multiple traces)

        FCy and FCySmooth are used similar to FCx, but could be recomputed
        from FCx if storage is more important than compute time.

        FCsta and FClta are not used in any way for processing, however they
        are useful for plotitng and visualizations.  Similar to FCy and
        FCySmooth, they could be computed from FCx if need be.

    So if the caches are cleared, which occurs when discontiguous traces are
    passed to calibration(), different snclq's are passed in or separate
    processes are instantiated, the caches will be cleared.  This means
    envelope normalization may not be accurate and if the calibration spans
    to earlier than the current trace it may not be identified. Also,
    filter initializations won't be accurate.  The filters to look out for
    are ones with long durations, such as the DC removal filter.  Short
    duration filters will adapt to the different initial conditions very
    quickly.

    ---------------------------------------------------------------------------

    :param trace:   an obspy trace object
    :param metric_store:    optional MetricStore
        (this is currently not used)
    :param debug:   flag, things are plotted and printed if set to True
    :param params: currently not used
    :return:  Some form of list of identified calibrations where each
    calibration is guaranteed to have:
        SNCLQ, start_time, end_time
    and may have additional statistics.  These details are specified in near
    the return statement.
    """

    def get_time_series_set(x):
        """Acquire all time-series' that will be used in detecting
        calibrations, primarily for initial trigger of candidate calibration
        detection.
        :param x: input time-series
        :return:
            x       - the input with DC removed (via LTI filter)
            xLP     - the input LP filtered to only contain DC
            sta     - short-term average
            lta     - long-term average
            y       - LTA / STA (that's right, the reciprocal of what you're
                       used to)
            ySmooth - filtered / smoothed version of 'y'
        """
        x = x.astype(dtype="float64")

        """remove DC"""

        if not _chelper.filters["DC_remover"].has_initial_condition():
            _chelper.filters["DC_remover"].zi = [np.mean(x)]

        xLP = _chelper.filters["DC_remover"].filter(x)
        x -= xLP

        """compute LTA/STA"""

        if not _chelper.filters["sta"].has_initial_condition():
            _chelper.filters["sta"].zi = [np.mean(abs(x))]
            _chelper.filters["lta"].zi = [np.mean(abs(x))]

        sta = _chelper.filters["sta"].filter(abs(x))
        lta = _chelper.filters["lta"].filter(abs(x))

        y = lta / sta

        """Smooth 'y' (maybe after detection check)"""

        if not _chelper.filters["ltasta_smoother"].has_initial_condition():
            _chelper.filters["ltasta_smoother"].zi = [np.mean(y)]

        ySmooth = _chelper.filters["ltasta_smoother"].filter(y)

        return (x, xLP, sta, lta, y, ySmooth)

    def get_candidate_indices(x):
        """indices where LTA/STA triggers above threshold"""

        # find values above threshold
        candidate_indices = [i for i in range(0, len(x)) if x[i] > _chelper.threshold]

        # remove trailing consecutive indices
        candidate_indices = [
            i for i in candidate_indices if i - 1 not in candidate_indices
        ]

        # remove 0 from indices if this is a continuation of a detection
        # from the previous frame or there is no previous frame
        if 0 in candidate_indices and (
            _chelper.FCySmooth.empty()
            or _chelper.FCySmooth.get(0)[-1] > _chelper.threshold
        ):
            del candidate_indices[0]

        return candidate_indices

    def initialize_detections(candidate_indices, ySmooth, trace):
        """Detections are managed as dictionaries where keys are added and
        deleted as necessary.  This function simply initializes this list
        of dictionaries."""
        cal_detections = []
        for idx in candidate_indices:
            temp = ySmooth[idx::]

            if min(temp) < _chelper.threshold:
                # ySmooth returns below threshold, find the snippet that is
                # above threshold (temp)
                idx2 = np.argmax([temp < _chelper.threshold]) - 1

                if idx2 == -1:
                    temp = temp[0:1]
                else:
                    temp = temp[0:idx2]

            # the initial end_time is at the trigger index
            end_time = trace.stats.starttime + idx / _chelper.fs

            cur_detection = {
                "idx": idx,
                "initial_cal_end_time": end_time,
                "max_ySmooth": max(temp),
            }
            cal_detections.append(cur_detection)

        return cal_detections

    def clear_failed_detections(cal_detections):
        for i in range(len(cal_detections) - 1, -1, -1):
            if "initial_cal_start_time" not in cal_detections[i]:
                del cal_detections[i]

    # flag that notifies the function to plot the data if debugging and a
    # candidate calibration fails some condition. Will never plot if 'debug'
    # is False
    _chelper.do_plot = False

    # update CalibrationHelper with new trace's information
    _chelper.update_state(trace)

    """generate all of the required time-series' """

    (x, xLP, sta, lta, y, ySmooth) = get_time_series_set(trace.data)

    """check for detections"""

    candidate_indices = get_candidate_indices(ySmooth)

    # if __calibration_demo__:
    #     import matplotlib.pyplot as plt
    #
    #     t = np.arange(len(x)) / trace.stats.sampling_rate
    #
    #     plt.subplot('211')
    #     plt.plot(t, x, color='k', linewidth=0.5,
    #              label='data with DC removed')
    #     plt.plot(t, xLP, color='r',
    #              label='DC that was removed')
    #     plt.plot(t, lta, color='y',
    #              label='Long Term Average (LTA)')
    #     plt.plot(t, sta, color='g',
    #              label='Short Term Average (STA)')
    #     plt.legend(loc=2)
    #     plt.xlim([t[0], t[-1]])
    #     plt.xlabel('t (seconds)')
    #     plt.ylabel('seismic data')
    #
    #     ax = plt.subplot('212');
    #     handles = []
    #     plt.plot(t, y, color='k', linewidth=0.5, linestyle='--',
    #              label='LTA / STA')
    #     plt.plot(t, ySmooth, color='k', linewidth=0.5,
    #              label='smoothed LTA / STA')
    #     plt.plot(t, _chelper.threshold * np.ones(len(t)),
    #              color='r', linestyle='--',
    #              label='threshold trigger level')
    #     plt.plot(t[candidate_indices], ySmooth[candidate_indices],
    #              linestyle='', marker='o', color='m',
    #              label='detection triggers')
    #
    #     if ax.get_ylim()[1] <= _chelper.threshold:
    #         ax.set_ylim([ax.get_ylim()[0], _chelper.threshold * 1.1])
    #     plt.legend(loc=2)
    #     plt.xlim([t[0], t[-1]])
    #     plt.xlabel('t (seconds)')
    #     plt.ylabel('LTA/STA data')
    #
    #     plt.show()

    """initialize the current detections"""

    cal_detections = initialize_detections(candidate_indices, ySmooth, trace)

    """Get start/end times for candidate detections"""

    if _chelper.FCzi.empty():
        # only do this once!
        _chelper.FCzi.add([np.mean(abs(x))])

    for cur_detection in cal_detections:
        idx = cur_detection["idx"]
        del cur_detection["idx"]

        """acquire accurate end time,
        (which is the transition where the LTA becomes > STA)"""

        y_reverse = y[idx:0:-1]

        # add previous frames until something under 1 is found
        for previous_y in _chelper.FCy:
            y_reverse = np.hstack((y_reverse, previous_y[::-1]))

            if min(y_reverse) < 1:
                break

        end_idx = np.argmax([y_reverse < 1])
        if end_idx == 0:
            # there is no transition point in all the data, void this
            # candidate index
            if debug:
                print("FAIL 1: finding accurate end_time")
                _chelper.do_plot = True
            continue
        end_idx = idx - end_idx

        issue_end_time = trace.stats.starttime + end_idx / _chelper.fs
        cur_detection["end_time"] = issue_end_time

        """find a start time"""

        absx_reverse = abs(x[max(end_idx, 0) : 0 : -1])

        # when looking for start, ignore indices below this
        search_start_idx = int(_chelper.minDetectionDuration * _chelper.fs)
        FCx_start_idx = 0

        # add frames until more than min detection duration is acquired
        for x_curframe in _chelper.FCx:
            FCx_start_idx += 1

            absx_rev_curframe = abs(x_curframe[::-1])

            # first, append frames until min search length is acquired
            absx_reverse = np.hstack((absx_reverse, absx_rev_curframe))

            if len(absx_reverse) > (search_start_idx + (end_idx < 0) * abs(end_idx)):
                break

        absx_reverse = absx_reverse[(end_idx < 0) * abs(end_idx) : :]

        if len(absx_reverse) <= search_start_idx:
            # not enough data to find a start, reject this candidate detection
            if debug:
                print("FAIL 2: not enough data for start_time")
                _chelper.do_plot = True
            continue

        # compute the envelope
        _chelper.filters["envelope"].zi = np.zeros(
            len(_chelper.filters["envelope"].B) - 1
        )
        env = _chelper.filters["envelope"].filter(absx_reverse)

        # normalize as a ratio to typical envelope level
        env /= _chelper.FCzi.get(_chelper.FCzi.size() - 1)[0]

        # add previous frames until something under 1 is found
        for cur_idx, x_curframe in enumerate(_chelper.FCx):
            if cur_idx < FCx_start_idx:
                # skip the frames that have already been added
                continue

            absx_rev_curframe = abs(x_curframe[::-1])

            # compute the envelope
            env_curframe = _chelper.filters["envelope"].filter(absx_rev_curframe)

            # normalize as a ratio to typical envelope level
            env_curframe /= _chelper.FCzi.get(_chelper.FCzi.size() - 1)[0]

            env = np.hstack((env, env_curframe))

            if min(env[search_start_idx::]) < 1:
                break

        start_idx = np.argmax([env[search_start_idx::] < 1]) + search_start_idx

        if start_idx == search_start_idx:
            # never found value below 1
            if not _chelper.FCx.full():
                # calibration extends to a data-not-available issue
                start_idx = len(env)
                cur_detection["do_recalibrate_start"] = False
            else:
                # no value below 1 was found within entire buffer,
                # assume that this was a false positive and discard
                if debug:
                    print("FAIL 3: envelope ratio never goes below 1")
                    _chelper.do_plot = True
                continue

        issue_start_time = issue_end_time - start_idx / _chelper.fs

        cur_detection["initial_cal_start_time"] = issue_start_time

        # prepare to collect envelope stats
        cur_detection["env"] = env[search_start_idx::]

    # delete any detections that failed the previous step
    clear_failed_detections(cal_detections)

    """recalibrate issue start times"""

    for cur_detection in cal_detections:
        issue_start_time = cur_detection["initial_cal_start_time"]
        issue_end_time = cur_detection["end_time"]

        if "do_recalibrate_start" in cur_detection:
            # issue extends to a data-not-available issue,
            # skip recalibration
            del cur_detection["do_recalibrate_start"]
            cur_detection["start_time"] = issue_start_time
            cur_detection["max_stalta"] = np.nan
        else:
            start_pad = _chelper.recal_tShort + _chelper.recal_tLong

            # get snippet of data from a bit before the calibration to the
            # end of the calibration
            temp_starttime = issue_start_time - start_pad
            temp_endtime = issue_end_time

            if trace.stats.starttime < temp_starttime:
                start_index = int(
                    (temp_starttime - trace.stats.starttime) * _chelper.fs
                )
                end_index = int((temp_endtime - trace.stats.starttime) * _chelper.fs)
                xtemp = trace.data[start_index:end_index]
            elif _chelper.fc_starttime < temp_starttime:
                start_index = int(
                    (temp_starttime - trace.stats.starttime) * _chelper.fs
                )
                end_index = int((temp_endtime - trace.stats.starttime) * _chelper.fs)

                xtemp = trace.data
                xtemp_start_index = 0
                for x_curframe in _chelper.FCx:
                    xtemp = np.hstack((x_curframe, xtemp))
                    xtemp_start_index -= len(x_curframe)

                    if xtemp_start_index <= start_index:
                        break

                start_index -= xtemp_start_index
                end_index -= xtemp_start_index
                try:
                    assert start_index >= 0 and start_index <= len(
                        xtemp
                    ), "not enough cached data"

                    xtemp = xtemp[start_index:end_index]
                except:
                    xtemp = None
            else:
                xtemp = None

            if xtemp is not None:
                # use the previously acquired snippent of time-series to
                # compute STA/LTA
                xtemp = xtemp.astype("float")

                xtemp -= np.mean(xtemp)

                xtemp = sta_lta(
                    xtemp,
                    _chelper.fs,
                    t_short=_chelper.recal_tShort,
                    t_long=_chelper.recal_tLong,
                )

                adjustment = np.argmax(xtemp) / _chelper.fs - 2 * _chelper.recal_tShort

                issue_start_time = issue_start_time + max(0, adjustment)

                cur_detection["start_time"] = issue_start_time
                cur_detection["max_stalta"] = max(xtemp)
            else:
                cur_detection["start_time"] = issue_start_time
                cur_detection["max_stalta"] = np.nan

        # select steady-state subset of envelope
        tempenv = cur_detection["env"]
        tempenv = tempenv[0 : int((issue_end_time - issue_start_time) * _chelper.fs)]
        tempenv = tempenv[0 : -np.argmax([tempenv[::-1] > np.median(tempenv)]) - 1]

        # save envelope statistics
        del cur_detection["env"]
        if len(tempenv) == 0:
            cur_detection["env_mean"] = np.nan
            cur_detection["env_median"] = np.nan
            cur_detection["env_var"] = np.nan
        else:
            cur_detection["env_mean"] = np.mean(tempenv)
            cur_detection["env_median"] = np.median(tempenv)
            cur_detection["env_var"] = np.var(tempenv)

    """store values to access in next frame"""

    # update time-series frame containers
    if _chelper.FCx.empty():
        _chelper.fc_starttime = trace.stats.starttime
    _chelper.FCx.add(x)
    _chelper.FCy.add(y)
    _chelper.FCsta.add(sta)
    _chelper.FClta.add(lta)
    _chelper.FCySmooth.add(ySmooth)
    _chelper.fc_endtime = trace.stats.endtime

    # update 'zi' frame container
    if len(cal_detections) == 0:
        # only update this if no detection occurred
        _chelper.filters["absxLP"].zi = _chelper.FCzi.get(0)
        _ = _chelper.filters["absxLP"].filter(abs(x))
        _chelper.FCzi.add(_chelper.filters["absxLP"].zi)
    else:
        # clear the previous frames that included a calibration
        # (this assumes the duration for each trace is constant)
        t = cal_detections[0]["start_time"]
        num_frames = int((trace.stats.starttime - t) / _chelper.bufferdur)
        zi = _chelper.FCzi.get(min(num_frames, _chelper.FCzi.size() - 1))
        _chelper.FCzi.clear()
        _chelper.FCzi.add(zi)

    """plot results if debugging or running demo"""
    # if _chelper.do_plot or __calibration_demo__:
    #     # add previous frames for plotting
    #     # num_extra_frames = 1
    #     num_extra_frames = 120
    #     num_extra_frames = min(num_extra_frames, _chelper.FCx.size() - 1)
    #     absx = abs(x)
    #     for i in range(1, num_extra_frames + 1):
    #         x = np.hstack((_chelper.FCx.get(i), x))
    #         absx = np.hstack((abs(_chelper.FCx.get(i)), absx))
    #         sta = np.hstack((_chelper.FCsta.get(i), sta))
    #         lta = np.hstack((_chelper.FClta.get(i), lta))
    #         y = np.hstack((_chelper.FCy.get(i), y))
    #         ySmooth = np.hstack((_chelper.FCySmooth.get(i), ySmooth))
    #
    #     import matplotlib.pyplot as plt
    #     from matplotlib import gridspec
    #     fig = plt.figure()
    #     gs = gridspec.GridSpec(2, 1)
    #     ax = fig.add_subplot(gs[0, 0])
    #     t = np.arange(0, len(x)) / _chelper.fs
    #     t -= num_extra_frames * _chelper.bufferdur
    #
    #     plt.plot(t, x, color='k', linewidth=.3, label='data with DC removed')
    #     plt.plot(t, _chelper.FCzi.get(_chelper.FCzi.size() - 1)[0]
    #              * np.ones(len(t)), color='r', linewidth=2.5,
    #              label='typical envelope level for normalization (DC value)')
    #     plt.plot(t, lta, color='y', linewidth=2,
    #              label='Long Term Average (LTA)')
    #     plt.plot(t, sta, color='g', linewidth=1.5,
    #              label='Short Term Average (STA)')
    #     plt.xlim([t[0], t[-1]])
    #     plt.ylabel('seismic data')
    #     plt.xlabel('t (seconds)')
    #     plt.legend(loc=2)
    #
    #     if len(cal_detections) > 0:
    #         cur_detection = cal_detections[-1]
    #         issue_end_time = cur_detection['end_time']
    #
    #         # plot envelope of last issue detection
    #         ax2 = ax.twinx()
    #
    #         t2 = np.arange(0, len(env)) / _chelper.fs
    #         t2 += issue_end_time - trace.stats.starttime - t2[-1]
    #         plt.plot(t2, 1 + 0 * env, color='b',
    #                  linestyle='--', linewidth=1,
    #                  label='envelope threshold (= 1)')
    #         plt.plot(t2, env[::-1], color='b',
    #                  linestyle='-', linewidth=3,
    #                  label='normalized envelope')
    #         plt.ylabel('envelope')
    #         plt.xlim([t[0], t[-1]])
    #         ax2.set_ylim([ax2.get_ylim()[1] * ax.get_ylim()[0] /
    #                       ax.get_ylim()[1],
    #                       ax2.get_ylim()[1]])
    #         plt.legend(loc=3)
    #
    #     if ax.get_ylim()[1] <= _chelper.threshold:
    #         ax.set_ylim([ax.get_ylim()[0], _chelper.threshold * 1.1])
    #     plt.xlabel('t (seconds)')
    #
    #     ax = fig.add_subplot(gs[1, 0])
    #     plt.plot(t, y, color='k', linestyle='--', linewidth=.8,
    #              label='LTA / STA')
    #     plt.plot(t, ySmooth, color='k', linewidth=1.5,
    #              label='smoothed LTA / STA')
    #     plt.plot(t, _chelper.threshold * np.ones(len(t)),
    #              color='r', linestyle='--', linewidth=3,
    #              label='threshold trigger level')
    #     plt.ylabel('LTA/STA data')
    #     plt.legend(loc=2)
    #
    #     # if len(issue_times) > 0:
    #     #     plt.plot(t[0:len(env)], env[::-1], color='y')
    #
    #     for cur_detection in cal_detections:
    #         ts = cur_detection['initial_cal_start_time']
    #         te = cur_detection['initial_cal_end_time']
    #
    #         ps = get_plot_specs()
    #         ps['linestyle'] = '--'
    #         ps['linewidth'] = 1
    #
    #         ps['color'] = 'g'
    #         issue_time = ts - trace.stats.starttime
    #         ax = fig.add_subplot(gs[0, 0])
    #         plot_vertical_line(ax, issue_time, ylims=[0, 10], ps=ps)
    #         ax = fig.add_subplot(gs[1, 0])
    #         plot_vertical_line(ax, issue_time, ylims=[0, 10], ps=ps)
    #
    #         ps['color'] = 'r'
    #         issue_time = te - trace.stats.starttime
    #         ax = fig.add_subplot(gs[0, 0])
    #         plot_vertical_line(ax, issue_time, ylims=[0, 10], ps=ps)
    #         ax = fig.add_subplot(gs[1, 0])
    #         plot_vertical_line(ax, issue_time, ylims=[0, 10], ps=ps)
    #
    #         ts = cur_detection['start_time']
    #         te = cur_detection['end_time']
    #
    #         ps = get_plot_specs()
    #         ps['linestyle'] = '-'
    #         ps['linewidth'] = 2
    #
    #         ps['color'] = 'g'
    #         issue_time = ts - trace.stats.starttime
    #         ax = fig.add_subplot(gs[0, 0])
    #         plot_vertical_line(ax, issue_time, ylims=[0, 10], ps=ps)
    #         ax = fig.add_subplot(gs[1, 0])
    #         plot_vertical_line(ax, issue_time, ylims=[0, 10], ps=ps)
    #
    #         ps['color'] = 'r'
    #         issue_time = te - trace.stats.starttime
    #         ax = fig.add_subplot(gs[0, 0])
    #         plot_vertical_line(ax, issue_time, ylims=[0, 10], ps=ps)
    #         ax = fig.add_subplot(gs[1, 0])
    #         plot_vertical_line(ax, issue_time, ylims=[0, 10], ps=ps)
    #
    #     plt.xlim([t[0], t[-1]])
    #     plt.ylim([ax.get_ylim()[0], max(ax.get_ylim()[1],
    #                                     1.1 * _chelper.threshold)])
    #     plt.legend(loc=2)
    #
    #     plt.show()

    """return final results"""

    return_detections = []
    for cur_detection in cal_detections:
        issue_start_time = cur_detection["start_time"]
        issue_end_time = cur_detection["end_time"]
        snclq = trace.get_id()

        """add additional statistics about the calibration"""

        # get the calibration segment of 'x'
        cal_x = []
        start_time = None
        for i in range(0, _chelper.FCx.size()):
            cal_x = np.hstack((_chelper.FCx.get(i), cal_x))
            if i == 0:
                start_time = trace.stats.starttime
            else:
                start_time -= len(_chelper.FCx.get(i)) / _chelper.fs
        start_idx = int((issue_start_time - start_time) * _chelper.fs)
        end_idx = int((issue_end_time - start_time) * _chelper.fs)
        # t = range(len(cal_x))
        # plt.plot(t, cal_x)
        # plt.plot(t[start_idx:end_idx], cal_x[start_idx:end_idx], color='r')
        # plt.show()
        cal_x = cal_x[start_idx:end_idx]

        # compute the STFT of the calibration
        window_length = 5  # seconds
        frame_size = window_length * _chelper.fs
        X, _, f_stft = stft(cal_x, fs=_chelper.fs, frame_size=frame_size)
        # get magnitudes of a subband
        flo = 0.5
        fhi = 20
        fidxlo = np.argmax([f_stft > flo])
        fidxhi = np.argmax([f_stft > fhi])
        if fidxhi == 0:
            fidxhi = len(f_stft)
        X = abs(X[:, fidxlo:fidxhi])
        # normalize
        X *= len(cal_x) ** 0.5 / np.sum(X ** 2) ** 0.5
        # get the variance of the STFT per frequency bin
        vars = np.var(X, axis=0)
        stft_var = np.mean(vars)

        """store results to detection"""

        rd = {  # start time (epoch time) of the identified calibration
            "cal_start_time": issue_start_time,
            # end time of the calibration
            "cal_end_time": issue_end_time,
            # initial detected start time before refining step
            "initial_cal_start_time": cur_detection["initial_cal_start_time"],
            # initial detected end time before refining step
            "initial_cal_end_time": cur_detection["initial_cal_end_time"],
            # maximum peak of initial detection trigger,
            # (larger implies more dramatic drop in energy)
            "max_ySmooth": cur_detection["max_ySmooth"],
            # maximum STA/LTA at start of calibration, similar to
            # max_ySmooth but for start, not end
            "max_stalta": cur_detection["max_stalta"],
            # average "energy" of the calibration
            # calibrations are likely to have a high env_mean
            "env_mean": cur_detection["env_mean"],
            # median "energy" of the calibration
            # calibrations are likely to have a high env_median
            "env_median": cur_detection["env_median"],
            # variance of the envelope
            # calibrations could have very high env_var but it's hard to say
            "env_var": cur_detection["env_var"],
            # mean variance of the STFT magnitude, if this is very high or
            # very low it is an indication of a calibration
            "stft_var": stft_var,
        }
        return_detections.append(rd)

    return return_detections


def calibration(trace, metric_store=None, debug=False, **params):
    """
    Refer to unscreened_calibration() for most details on this function.
    This function wraps unscreened_calibration and discards false alarms.
    At the moment it uses pre-trained machine learning objects however this
    can be implemented in different ways in the future.  The keys in a
    detection that are likely of interest to a user include:
      'snclq', 'start_time' and 'end_time'
    Beyond that, these entries are strictly for screening false alarms.
    """

    detections = unscreened_calibration(trace, metric_store, debug, **params)

    for i in range(len(detections) - 1, -1, -1):
        if not _chelper.iscal(detection=detections[i]):
            del detections[i]
        # Remove detections that are triggering on beginning of trace
        try:
            if detections[i]["cal_start_time"] == trace.stats.starttime:
                del detections[i]
        except IndexError:
            continue
    # Remove identical detections -- TODO this should really be fixed in the threshold triggering
    if len(detections) != 0:
        iden_detections = [d["cal_start_time"] for d in detections]
        dups = [
            idx
            for idx, item in enumerate(iden_detections)
            if item in iden_detections[:idx]
        ]
        for i in range(len(dups)):
            detections[dups[i]] = np.nan
        detections = [x for x in detections if str(x) != "nan"]

    # converting time to isoformat
    for i in range(len(detections)):
        for key, value in detections[i].items():
            if "time" in key:
                detections[i][key] = value.isoformat()

    # If no detections are returned, create the following list dictionary
    out = {
        "snclq": trace.get_id(),
        "start_time": trace.stats.starttime.isoformat(),
        "end_time": trace.stats.endtime.isoformat(),
        "num_cals_detected": len(detections),
        "metric_name": "calibrationMetric",
        "detections": detections,
    }

    return out


def calibrationMetric(st, metric_store=None, debug=False, database=None, **params):
    """
    Wrapper function for calibration. Enable stream processing
    :param st: Stream object

    :return: List of calibration detections
    """

    out = []
    for i in range(len(st)):
        tr = st[i]
        d = calibration(tr, metric_store, debug, **params)
        out.append(d)

    if database is not None or isinstance(st, Database):
        database.insert_metric(out)

    return out


# if __name__ == '__main__':
#
#     import os
#
#     def print_results(detections):
#         for detection in detections:
#             print('snclq            : %s' % detection['snclq'])
#             print('issue_start_time : %.2f' % detection['cal_start_time'])
#             print('issue_end_time   : %.2f' % detection['cal_end_time'])
#             for key in ['initial_cal_start_time',
#                         'cal_start_time',
#                         'cal_end_time',
#                         'initial_cal_end_time',
#                         'max_ySmooth',
#                         'max_stalta',
#                         'env_mean',
#                         'env_median',
#                         'env_var',
#                         'stft_var']:
#                 print('%20s  :  %s' % (key, str(detection[key])))
#
#
#     __calibration_demo__ = False
#
#     # test single trace
#     data = os.path.dirname(os.path.dirname(__file__)) + ''
#     # st = obspy.read(data)
#     # tr = st[0]
#     data = '/Users/prkay/Workspace/pycheron/pycheron/metrics/test_data/ZNPU_HHN_006.mseed'
#
#     detections = calibrationMetric(obspy.read(data))
#     print(detections)
#     print_results(detections)
