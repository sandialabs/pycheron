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

import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path + "/..")
import warnings
from obspy import UTCDateTime
import pandas as pd
import numpy as np
from calendar import monthrange

__all__ = ["MetricStore"]


class MetricStore(object):
    """
    This class contains average stats for each unique snclq.  The stats are
    possibly metric values but may also be statistics important to the
    metric function but hidden from the user. metric values may be a number
    or an array-like sequence of numbers. This can vary but it needs to be
    consistent for any specific SNCLQ/metric_name combination.

    Data is updated as 1-day values, which are stored in local memory. The
    single-day data is gathered and once a full month of data is in local
    memory, the data is averaged and stored to non-volatile memory.  This is
    designed to be used in a streaming manner, as can be seen in the example in
    __main__.  If a discontiguous date or SNCLQ is added to the metric
    store, the previous per-day data is cleared.

    The actual files are handles as follows:
    -new file per SNCL
    -new DataFrame per metric
    -new row per year and month
    """

    def __init__(self, id=0):
        """
        :param id: Unique ID for a store (integer)
        Each time a new id is provided a new subdirectory is created where
        that store's data is contained.  A new store may be created if,
        for example, a new version of the code is developed it is desired
        that the old data be left untouched, or multiple datasets are are
        being processed and it is desired that their data be kept separate.
        """
        self._id = id
        # self._storedir = dir_path + '/streaming/data/%03d/' % self._id
        self._storedir = dir_path + "/data/%03d/" % self._id
        if not os.path.exists(self._storedir):
            os.makedirs(self._storedir)
        self._day_tracker = self._DayTracker()

    def metric_names(self, snclq):
        """
        Returns the names of each metric that has already been added to the
        metric store for this SNCLQ
        :param snclq:
        :return: list of strings
        """
        storepath = self._storepath(snclq)
        if not os.path.exists(storepath):
            return []

        with pd.HDFStore(storepath, mode="r") as store:
            names = list(store.keys())

        # remove the '/' prepended by HDF5 store
        for i in range(len(names)):
            names[i] = names[i].lstrip("/")

        return names

    def put(
        self, snclq, metric_name, metric_value, date, overwrite=False, reshape=False
    ):
        """
        put new metric value into local memory, then if the date is the last
        day in the month and the rest of the month is in local memory,
        write the full month's data to the store.  The data written to the
        store for a month is the average of the values for each day.

        :param snclq:
        :param metric_name: a unique metric name;
            not necessarily something that the user would need!
        :param metric_value: new value (must be a number or array-like)
        :param overwrite: boolean; if False, the non-volatile store will not be
            modified
        :return: True on successful write to store

        NOTE: it is important that the metric_value is a consistent
        dimension for each unique metric, or averaging will throw an error.
        """

        if self._day_tracker.add(snclq, metric_name, metric_value, date):
            self._update(snclq, metric_name, metric_value, overwrite, reshape)

    def get_days(self, snclq, metric_name):
        """
        Get the per-day data currently stored in volatile memory. Each time
        a discontinuous date is added to the store, the previously stored
        days will be cleared.
        :param snclq:
        :param metric_name: unique name for a metric
        :return: tuple of lists: (values, dates)
         "values" is a list of numbers or arrays, whatever was provided as
            metric_value
         "dates" is a list of datetime objects that are always quantized to
            day, i.e. hour=0, minute=0, etc.
        """
        return self._day_tracker.get_stored_data(snclq, metric_name)

    def get(self, snclq, metric_name, year_month):
        """
          Get the value of a metric for a specified snclq and year_month
        :param snclq:
        :param metric_name: string name of metric
        :param year_month: tuple of two integers: (year, month)
            OR any object that contains "year" and "month" attributes
        :return: metric value or None if entry is missing
            metric value is an array of numbers or a number
        """
        data = None
        storepath = self._storepath(snclq)

        if os.path.exists(storepath):
            with pd.HDFStore(storepath, mode="r") as store:
                if metric_name in store:
                    ym_index = self._year_month_index(year_month)

                    try:
                        df = store[metric_name]
                    except TypeError:
                        del store[metric_name]
                        assert False, "data corrupted!"
                    except KeyError:
                        df = None

                    if df is not None and ym_index in df.index:
                        data = df.loc[ym_index].values

        return data

    """private methods"""

    class _DayTracker(object):
        """Object that manages the per-day data that is stored in volatile
        memory before the average for a whole month is written to the store
        in non-volatile memory."""

        def __init__(self, max_size=31):
            self._flush("", "")
            self._max_size = max_size

        def add(self, snclq, metric_name, metric_value, date):
            """
              Add a new entry to the tracker
            :param snclq:
            :param metric_name:
            :param metric_value:
            :param date: UTCDateTime or something that can be converted to one
            :return: True or False
                True implies a new complete month of data is available,
                i.e. the last day of the month was just added and the rest
                of the month is in the cache
            """
            date = self._quantized_day(UTCDateTime(date))

            if not self._is_same(snclq, metric_name) or (
                len(self.dates) > 0 and date - self.dates[-1] != 24 * 3600
            ):
                """
                either:
                   1) the snclq or metric_name is inconsistent
                   2) the new date is not the day after the last stored date
                therefore, flush the cached days
                """
                self._flush(snclq, metric_name)

            # check that metric_value dimensions are consistent
            if len(self.values) > 0:
                assert len(self.values[0]) == len(metric_value), (
                    "inconsistent metric_value length\nexpected: "
                    "%d\nactual: %d\n" % (len(self.values[0]), len(metric_value))
                )

            # add new data entry
            self.values.append(metric_value)
            self.dates.append(date)

            # remove oldest data if past cache size
            if len(self.dates) > self._max_size:
                self.values = self.values[-self._max_size : :]
                self.dates = self.dates[-self._max_size : :]

            # new month of data ready?
            return (date + 24 * 3600).month != date.month and len(
                self.dates
            ) >= date.day

        def get_stored_data(self, snclq, metric_name):
            if self._is_same(snclq, metric_name):
                return (self.values, self.dates)
            else:
                return ([], [])

        def get_year_month(self):
            if len(self.dates) == 0:
                return None
            else:
                return (self.dates[-1].year, self.dates[-1].month)

        def get_year_month_data(self):
            if self.get_year_month() is None:
                return None, None
            else:
                (year, month) = self.get_year_month()

            data = [
                self.values[i]
                for i, date in enumerate(self.dates)
                if date.year == year and date.month == month
            ]

            if len(data) != monthrange(year, month)[1]:
                # not a complete month of data!
                return None, None

            return data, (year, month)

        def _is_same(self, snclq, metric_name):
            return snclq == self.snclq and metric_name == self.metric_name

        def _flush(self, snclq, metric_name):
            self.snclq = snclq
            self.metric_name = metric_name
            self.values = []
            self.dates = []

        @staticmethod
        def _quantized_day(datetime):
            return datetime - datetime.timestamp % (24 * 3600)

    def _storepath(self, snclq):
        """returns the path to the SNCLQ-specific store"""
        return self._storedir + snclq + ".h5"

    @staticmethod
    def _year_month_index(year_month):
        """
        Converts a year_month to a index string and also checks for valid input
        :param year_month: either:
            1) an object with 'year' and 'month' attributes
            2) a tuple of the format (year, month), both as integers
        :return: index string for accessing into an HDF5 store
        """
        if hasattr(year_month, "year") and hasattr(year_month, "month"):
            year = year_month.year
            month = year_month.month
        else:
            assert len(year_month) == 2, (
                "invalid input year_month"
                "\nmust be datetime-like object or "
                "tuple containing (year, month)"
            )
            (year, month) = year_month

        # check that year and month values are valid
        _ = UTCDateTime(year, month, 1)

        return "%4d-%2d" % (year, month)

    def _update(self, snclq, metric_name, metric_value, overwrite=False, reshape=False):
        """Write or re-write the entry for the specified snclq and
        metric_name and year_month. Uses the current day-data in the day
        tracker.  This function is called every time the new entry in the
        day tracker completes the month"""

        # data ready to be written to store
        new_data, year_month = self._day_tracker.get_year_month_data()

        if not overwrite and self.get(snclq, metric_name, year_month) is not None:
            # data already written, don't overwrite
            return

        # average the data over all days in the month
        avg_data = np.mean(np.array(new_data), axis=0)

        ym_index = self._year_month_index(year_month)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pd.HDFStore(self._storepath(snclq), mode="a") as store:
                if (
                    metric_name not in store
                    # or... the new metric_value and stored metric_value
                    # have inconsistent lengths
                    or (reshape and len(avg_data) != store[metric_name].shape[1])
                ):
                    # first entry for this metric
                    df = pd.DataFrame([avg_data], index=[ym_index])
                    store[metric_name] = df
                else:
                    try:
                        # overwrite or append index depending on if it is already
                        # in the data frame
                        assert len(avg_data) == store[metric_name].shape[1], (
                            "new metric_value and stored metric_value have "
                            "inconsistent lengths but reshape=False! Either "
                            "call with reshape=True or make sure the "
                            "dimensions are consistent. WARNING: reshape "
                            "will clear any stored values for the current "
                            "snclq/metric_name"
                            "\nLength in store: %d"
                            "\nLength of new:   %d"
                            % (store[metric_name].shape[1], len(avg_data))
                        )
                        store[metric_name].loc[ym_index] = avg_data
                    except TypeError:
                        del store[metric_name]
                        assert False, "data corrupted!"


class MultiSnclqMetricStore(object):
    """
    This is a copy-paste of the MetricStore function that Kale and I are
    going to modify in the middle of a meeting. NOTE: none of this has been
    tested as of 6/6/17 but it should "theoretically" work fine.

    Functionality:
    This class operates the exact same way as MetricStore, and in this sense
    docstrings for MetricStore may be used to help understand its functionality
    WITH THE EXCEPTION that instead of having a single day_tracker that
    caches day-data, there will be a new day_tracker for each unique SNCLQ.
    This means that it will take up more volatile memory, but that adding
    new data from one snclq won't imply clearing the cache of day-data for
    the other snclq's
    """

    def __init__(self, id=0):
        """
        :param id: Unique ID for a store (integer)
        Each time a new id is provided a new subdirectory is created where
        that store's data is contained.  A new store may be created if,
        for example, a new version of the code is developed it is desired
        that the old data be left untouched, or multiple datasets are are
        being processed and it is desired that their data be kept separate.
        """
        dir_path = os.path.dirname(os.path.realpath(__file__))

        self._id = id
        # self._storedir = dir_path + 'streaming/data/%03d/' % self._id
        self._storedir = dir_path + "/data/%03d/" % self._id
        print(self._storedir)

        if not os.path.exists(self._storedir):
            os.makedirs(self._storedir)

        self._day_trackers = {}

    def metric_names(self, snclq):
        """
        Returns the names of each metric that has already been added to the
        metric store for this SNCLQ
        :param snclq:
        :return: list of strings
        """
        storepath = self._storepath(snclq)
        if not os.path.exists(storepath):
            return []

        with pd.HDFStore(storepath, mode="r") as store:
            names = list(store.keys())

        # remove the '/' prepended by HDF5 store
        for i in range(len(names)):
            names[i] = names[i].lstrip("/")

        return names

    def put(
        self, snclq, metric_name, metric_value, date, overwrite=False, reshape=False
    ):
        """
        put new metric value into local memory, then if the date is the last
        day in the month and the rest of the month is in local memory,
        write the full month's data to the store.  The data written to the
        store for a month is the average of the values for each day.

        :param snclq:
        :param metric_name: a unique metric name;
            not necessarily something that the user would need!
        :param metric_value: new value (must be a number or array-like)
        :param overwrite: boolean; if False, the non-volatile store will not be
            modified
        :return: True on successful write to store

        NOTE: it is important that the metric_value is a consistent
        dimension for each unique metric, or averaging will throw an error.
        """
        day_tracker = self._day_tracker(snclq)

        if day_tracker.add(snclq, metric_name, metric_value, date):
            self._update(snclq, metric_name, metric_value, overwrite, reshape)

    def get_days(self, snclq, metric_name):
        """
        Get the per-day data currently stored in volatile memory. Each time
        a discontinuous date is added to the store, the previously stored
        days will be cleared.
        :param snclq:
        :param metric_name: unique name for a metric
        :return: tuple of lists: (values, dates)
         "values" is a list of numbers or arrays, whatever was provided as
            metric_value
         "dates" is a list of datetime objects that are always quantized to
            day, i.e. hour=0, minute=0, etc.
        """
        day_tracker = self._day_tracker(snclq)
        return day_tracker.get_stored_data(snclq, metric_name)

    def get(self, snclq, metric_name, year_month):
        """
          Get the value of a metric for a specified snclq and year_month
        :param snclq:
        :param metric_name: string name of metric
        :param year_month: tuple of two integers: (year, month)
            OR any object that contains "year" and "month" attributes
        :return: metric value or None if entry is missing
            metric value is an array of numbers or a number
        """
        data = None
        storepath = self._storepath(snclq)

        if os.path.exists(storepath):
            with pd.HDFStore(storepath, mode="r") as store:
                if metric_name in store:
                    ym_index = self._year_month_index(year_month)

                    try:
                        df = store[metric_name]
                    except TypeError:
                        del store[metric_name]
                        assert False, "data corrupted!"
                    except KeyError:
                        df = None

                    if df is not None and ym_index in df.index:
                        data = df.loc[ym_index].values

        return data

    """private methods"""

    class _DayTracker(object):
        """Object that manages the per-day data that is stored in volatile
        memory before the average for a whole month is written to the store
        in non-volatile memory."""

        def __init__(self, max_size=31):
            self._flush("", "")
            self._max_size = max_size

        def add(self, snclq, metric_name, metric_value, date):
            """
              Add a new entry to the tracker
            :param snclq:
            :param metric_name:
            :param metric_value:
            :param date: UTCDateTime or something that can be converted to one
            :return: True or False
                True implies a new complete month of data is available,
                i.e. the last day of the month was just added and the rest
                of the month is in the cache
            """
            date = self._quantized_day(UTCDateTime(date))

            if not self._is_same(snclq, metric_name) or (
                len(self.dates) > 0 and date - self.dates[-1] != 24 * 3600
            ):
                """
                either:
                   1) the snclq or metric_name is inconsistent
                   2) the new date is not the day after the last stored date
                therefore, flush the cached days
                """
                self._flush(snclq, metric_name)

            # check that metric_value dimensions are consistent
            if len(self.values) > 0:
                assert len(self.values[0]) == len(metric_value), (
                    "inconsistent metric_value length\nexpected: "
                    "%d\nactual: %d\n" % (len(self.values[0]), len(metric_value))
                )

            # add new data entry
            self.values.append(metric_value)
            self.dates.append(date)

            # remove oldest data if past cache size
            if len(self.dates) > self._max_size:
                self.values = self.values[-self._max_size : :]
                self.dates = self.dates[-self._max_size : :]

            # new month of data ready?
            return (date + 24 * 3600).month != date.month and len(
                self.dates
            ) >= date.day

        def get_stored_data(self, snclq, metric_name):
            if self._is_same(snclq, metric_name):
                return (self.values, self.dates)
            else:
                return ([], [])

        def get_year_month(self):
            if len(self.dates) == 0:
                return None
            else:
                return (self.dates[-1].year, self.dates[-1].month)

        def get_year_month_data(self):
            if self.get_year_month() is None:
                return None, None
            else:
                (year, month) = self.get_year_month()

            data = [
                self.values[i]
                for i, date in enumerate(self.dates)
                if date.year == year and date.month == month
            ]

            if len(data) != monthrange(year, month)[1]:
                # not a complete month of data!
                return None, None

            return data, (year, month)

        def _is_same(self, snclq, metric_name):
            return snclq == self.snclq and metric_name == self.metric_name

        def _flush(self, snclq, metric_name):
            self.snclq = snclq
            self.metric_name = metric_name
            self.values = []
            self.dates = []

        @staticmethod
        def _quantized_day(datetime):
            return datetime - datetime.timestamp % (24 * 3600)

    def _storepath(self, snclq):
        """returns the path to the SNCLQ-specific store"""
        return self._storedir + snclq + ".h5"

    @staticmethod
    def _year_month_index(year_month):
        """
        Converts a year_month to a index string and also checks for valid input
        :param year_month: either:
            1) an object with 'year' and 'month' attributes
            2) a tuple of the format (year, month), both as integers
        :return: index string for accessing into an HDF5 store
        """
        if hasattr(year_month, "year") and hasattr(year_month, "month"):
            year = year_month.year
            month = year_month.month
        else:
            assert len(year_month) == 2, (
                "invalid input year_month"
                "\nmust be datetime-like object or "
                "tuple containing (year, month)"
            )
            (year, month) = year_month

        # check that year and month values are valid
        _ = UTCDateTime(year, month, 1)

        return "%4d-%2d" % (year, month)

    def _update(self, snclq, metric_name, metric_value, overwrite=False, reshape=False):
        """Write or re-write the entry for the specified snclq and
        metric_name and year_month. Uses the current day-data in the day
        tracker.  This function is called every time the new entry in the
        day tracker completes the month"""

        # data ready to be written to store
        day_tracker = self._day_tracker(snclq)
        new_data, year_month = day_tracker.get_year_month_data()

        if not overwrite and self.get(snclq, metric_name, year_month) is not None:
            # data already written, don't overwrite
            return

        # average the data over all days in the month
        avg_data = np.nanmean(np.array(new_data), axis=0)

        ym_index = self._year_month_index(year_month)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pd.HDFStore(self._storepath(snclq), mode="a") as store:
                if (
                    metric_name not in store
                    # or... the new metric_value and stored metric_value
                    # have inconsistent lengths
                    or (reshape and len(avg_data) != store[metric_name].shape[1])
                ):
                    # first entry for this metric
                    df = pd.DataFrame([avg_data], index=[ym_index])
                    store[metric_name] = df
                else:
                    try:
                        # overwrite or append index depending on if it is already
                        # in the data frame
                        assert len(avg_data) == store[metric_name].shape[1], (
                            "new metric_value and stored metric_value have "
                            "inconsistent lengths but reshape=False! Either "
                            "call with reshape=True or make sure the "
                            "dimensions are consistent. WARNING: reshape "
                            "will clear any stored values for the current "
                            "snclq/metric_name"
                            "\nLength in store: %d"
                            "\nLength of new:   %d"
                            % (store[metric_name].shape[1], len(avg_data))
                        )
                        store[metric_name].loc[ym_index] = avg_data
                    except TypeError:
                        del store[metric_name]
                        assert False, "data corrupted!"

    def _day_tracker(self, snclq):
        """
          Acquire (and create if necessary) a _DayTracker object for the
          provided snclq. This way day data can be cached for each snclq
          without worrying about flushing data for the other snclq's.
        :param snclq:
        :return: a _DayTracker object
        """
        if snclq in self._day_trackers:
            # this is a previously input SNCLQ, use the already-generated
            # day tracker
            day_tracker = self._day_trackers[snclq]
        else:
            # this is not a previously input SNCLQ, create a new day tracker
            #  for it
            day_tracker = self._DayTracker()
            self._day_trackers[snclq] = day_tracker

        return day_tracker


# if __name__ == '__main__':
#     # initialize the metric store
#     ms1 = MetricStore(id=0)
#     ms0 = MultiSnclqMetricStore(id=0)
#     ms_list = [ms0, ms1]
#     for ms in ms_list:
#         # specify whether or not we should overwrite the metric store data
#         # permanently
#         overwrite = True
#
#         # specify a unique snclq
#         snclq = 'LPAZ.BHZ'
#         print(ms.metric_names(snclq))
#
#         # specify a unique metric_name
#         metric_name = 'calibration.env_var'
#
#         # get that metric value for a specific year and month
#         year_month = (2012, 9)
#         print(repr(ms.get(snclq, metric_name, year_month)))
#
#         for number_of_days in range(60):
#             # move a day forward
#             date = UTCDateTime(2012, 8, 28) + (24 * 3600 * number_of_days)
#
#             # create arbitrary data
#             metric_value = np.random.randn(3)
#             # if you try this instead, the inconsistent dimension will throw an
#             # error
#             # metric_value = np.random.randn(5)
#
#             ms.put(snclq, metric_name, metric_value, date,
#                    overwrite=overwrite, reshape=False)
#
#             if number_of_days % 10 == 0:
#                 print('\nData in local memory:')
#                 values, dates = ms.get_days(snclq, metric_name)
#                 for i in range(len(values)):
#                     print('%2d/%2d/%4d  -  %s' % (dates[i].month, dates[i].day,
#                                                   dates[i].year, str(values[i])))
#
#         print(repr(ms.get(snclq, metric_name, year_month)))
