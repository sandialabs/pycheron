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

__all__ = ["evaluate", "evaluate_stream", "evaluate_file", "evaluate_webservice"]

import os
from obspy import read
from obspy.clients.fdsn import Client
import joblib
import numpy as np
import os
import sklearn
from pycheron.sigpro.qcStatisticsML.feature_generator import FeatureGenerator
from pycheron.db.sqllite_db import Database


def evaluate(
    trace,
    model_path="/sigpro/qcStatisticsML/models/1621881059_RF-Raw-QC_Balanced.joblib.pkl",
    window_size=60,
    stride_length=50,
):
    """

    Rolling window function to evaluate whether observed Obspy trace object contains signal or an artifact using a
    Random Forest ML model

    :param tr: obspy trace object
    :type tr: obspy.core trace object
    :param: model_path: training model file path location
    :type model_path: str
    :param window_size: rolling window size in seconds (default = 60s)
    :type window_size: int
    :param stride_length: length in seconds to increment rolling window (default = 50s)
    :type stride_length: int

    returns: d; dictionary containing all the artifacts found for each trace object. The following keys exist:

            * snclq (`str`): station network channel location quality code information for specfic trace
            * start_time (`str`): start time of trace
            * end_time (`str`): end time of trace
            * metric_name (`str`): qcMLMetric string
            * qc_ml_results (`list`): list of dictionaries containing the following key/values
                * window_start: (`UTCDateTime`) - start time of rolling window
                * window_end: (`UTCDateTime`) - end time of rolling window
                * Features: (`numpy array`) - array of feature values: drop out fraction, distinct values ratio, packet time
                                            bandwidth product, frequency sigma, and discontinuity max value.
                                            * dropout fraction is the number of discrete intervals in the trace with N
                                                or more consecutive samples having the same value
                                            * distinct values is the number of distinct values divided by the total
                                                number of values (this is inversely related to the quantization error)
                                            * tbp is packet time bandwidth product
                                            * frequency sigma - standard deviation of the frequency
                                            * discontinuity max value - maximum value of discontinuities
                * prediction": (`numpy.array`) - array of a list with a string. Since only concerned with artifacts, these
                                                will all have the following value ['artifact']

    :rtype: dict

    * Code base originally written by Steven Magana-Zook from Lawrence Livermore National Laboratories
      Augmented and adapted for use within Pycheron by Pycheron team

    ** Resources for ML/features selected Method **:

    #. Dodge, D.A., and D.B. Harris (2016). Local-Scale Test of Dynamic Correlation Processors: Implications for
       Correlation-Based Seismic Pipelines, Bull Seis. Am. 106(2), pp. 435-452.

    #. Dodge, D.A., and W.R. Walter (2015). Initial Global Seismic Cross-Correlation Results: Implications for Empirical
       Signal Detectors, Bull. Seis. Am. 105(1), pp.240-256.

    """

    # Get current working directory to update model path to full path.
    # Do we need the models to be in Pycheron? If so, keep model_path.
    # If not, possibly get rid of this in the future and have user
    # paste full path of file
    model_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + model_path

    # Set up list to append to when find artifacts. Add basic trace information
    d = {
            "snclq": trace.get_id(),
            "start_time": trace.stats.starttime,
            "end_time": trace.stats.endtime,
            "metric_name": "qcMLMetric",
            "dropout_fraction": [],
            "distinct_values_ratio": [],
            "packet_time_bandwidth_product": [],
            "frequency_sigma": [],
            "discontinuity_max_value": [],
            "artifacts": 0,
            "qc_ml_results": []
        }

    # Load pickle ML model file
    clf_rf = joblib.load(model_path)

    # Call feature generator class to calculate drop out fraction, distinct values, tbp, frequency sigma and
    # disc_max_value features
    f_generator = FeatureGenerator()

    # Set window start time to the start time of the trace and then iterate through the trace object
    window_start_time = trace.stats.starttime
    # While the window start time + size of the window is less than the endtime of the trace, evaluate each window
    # segment against ML model and determine whether the result was a signal or an artifact
    while window_start_time + window_size <= trace.stats.endtime:
        # Get end time of window segment, then cut the trace based on start/end time to get isolated evaluation window
        window_end_time = window_start_time + window_size
        eval_segment_window = trace.slice(window_start_time, window_end_time)

        # Run this window segment through the ML model and make predictions
        features = f_generator.generate(eval_segment_window)
        features_one_instance = np.array(features).reshape(1, -1)
        res = clf_rf.predict(features_one_instance)

        # Capture results if flagged as artifact
        # Ensuring use of native python types so results are JSON serializable during db ingestion
        if res[0] == "artifact":
            data = {
                "window_start": window_start_time.isoformat(),
                "window_end": window_end_time.isoformat(),
                "Features": features_one_instance.tolist(),
                "prediction": res.tolist(),
            }
            # Append result to overall list
            d["dropout_fraction"].append(features_one_instance[0][0])
            d["distinct_values_ratio"].append(features_one_instance[0][1])
            d["packet_time_bandwidth_product"].append(features_one_instance[0][2])
            d["frequency_sigma"].append(features_one_instance[0][3])
            d["discontinuity_max_value"].append(features_one_instance[0][4])
            d["artifacts"] += 1
            d["qc_ml_results"].append(data)

        # increment window using stride length
        window_start_time = window_start_time + stride_length
    pass

    # If no artifacts identified fill in prediction key with string 'No artifacts found'
    if len(d) == 1 and "prediction" not in d[0]:
        d[0].update({"prediction": "No artifacts found"})

    return d


def evaluate_stream(
    st,
    model_path="/sigpro/qcStatisticsML/models/1715965925_RF-Raw-QC_Balanced.joblib.pkl",
    window_size=60,
    stride_length=50,
    database_config=None,
):
    """
    Function to loop through each traces in stream object and evaluate whether it contains signal or an artifact using
    a random forest ML model

    :param st: obspy stream object
    :type st: obspy.core stream object
    :param: model_path: training model file path location. Should be relative location as evaluate function will look for
                        full directory path. Default example(
                        '/sigpro/qcStatisticsML/models/1621881059_RF-Raw-QC_Balanced.joblib.pkl')
    :type model_path: str
    :param window_size: rolling window size in seconds
    :type window_size: int
    :param stride_length: length in seconds to increment rolling window
    :type stride_length: int
    :param database: Pycheron database to ingest metrics into
    :type database: pycheron.db.sqllite_db.Database

    returns: artifact_list; list of dictionaries containing all the artifacts found for each trace object within
                            respective stream object. The following keys exist:

            * snclq (`str`): station network channel location quality code information for specfic trace
            * start_time (`str`): start time of trace
            * end_time (`str`): end time of trace
            * metric_name (`str`): qcMLMetric string
            * qc_ml_results (`list`): list of dictionaries containing the following key/values
                * window_start: (`UTCDateTime`) - start time of rolling window
                * window_end: (`UTCDateTime`) - end time of rolling window
                * Features: (`numpy array`) - array of feature values: drop out fraction, distinct values ratio, packet time
                                            bandwidth product, frequency sigma, and discontinuity max value.
                                            * dropout fraction is the number of discrete intervals in the trace with N
                                                or more consecutive samples having the same value
                                            * distinct values is the number of distinct values divided by the total
                                                number of values (this is inversely related to the quantization error)
                                            * tbp is packet time bandwidth product
                                            * frequency sigma - standard deviation of the frequency
                                            * discontinuity max value - maximum value of discontinuities
                * prediction": (`numpy.array`) - array of a list with a string. Since only concerned with artifacts, these
                                                will all have the following value ['artifact']
            * dropout_fraction (`list`): series of results from the feature metric, starts at first window/ends at last
            * distinct_values_ratio (`list`): series of results from the feature metric, starts at first window/ends at last
            * packet_time_bandwidth_product (`list`): series of results from the feature metric, starts at first window/ends at last
            * frequency_sigma (`list`): series of results from the feature metric, starts at first window/ends at last
            * discontinuity_max_value (`list`): series of results from the feature metric, starts at first window/ends at last
            * artifacts (`list`): series of results from the feature metric, starts at first window/ends at last

    :rtype: list of dictionaries
    """

    # Set up list to append output to for each trace
    artifact_list = []
    # Loop through traces in stream object and evaluate each trace for artifacts. Append to artifact_list to build up
    # list of artifacts for stream object
    for trace in st:
        artifact_list.append(
            evaluate(
                trace,
                model_path=model_path,
                window_size=window_size,
                stride_length=stride_length,
            )
        )

    if database_config is not None:
        database = Database(**database_config)
        database.insert_metric(artifact_list)

    return artifact_list


def evaluate_file(
    file_path,
    model_path="/sigpro/qcStatisticsML/models/1621881059_RF-Raw-QC_Balanced.joblib.pkl",
    window_size=60,
    stride_length=50,
    database_config=None,
):
    """
    Function to evaluate a user-provided file, convert it to a stream object, and loop through each trace in the
    stream object (using evaluate_stream function) to evaluate whether it contains signal or an artifact using a Random
    Forest ML model (using evaluate function)

    :param file_path: input file path location
    :type file_path: str
    :param: model_path: training model file path location. Should be relative location as evaluate function will look for
                        full directory path. Default example(
                        '/sigpro/qcStatisticsML/models/1621881059_RF-Raw-QC_Balanced.joblib.pkl')
    :type model_path: str
    :param window_size: rolling window size in seconds
    :type window_size: int
    :param stride_length: length in seconds to increment rolling window
    :type stride_length: int
    :param database_config: dictionary containing the necessary parameters to create
                            a pycheron Database object. 
                            These include "db_name", "session_name", "overwrite", "manual", "wfdb_conn"
    :type database_config: dict

     returns: d; list of dictionaries containing all the artifacts found for each trace object within
                            respective stream object. The following keys exist:

            * snclq (`str`): station network channel location quality code information for specfic trace
            * start_time (`str`): start time of trace
            * end_time (`str`): end time of trace
            * metric_name (`str`): qcMLMetric string
            * qc_ml_results (`list`): list of dictionaries containing the following key/values
                * window_start: (`UTCDateTime`) - start time of rolling window
                * window_end: (`UTCDateTime`) - end time of rolling window
                * Features: (`numpy array`) - array of feature values: drop out fraction, distinct values ratio, packet time
                                            bandwidth product, frequency sigma, and discontinuity max value.
                                            * dropout fraction is the number of discrete intervals in the trace with N
                                                or more consecutive samples having the same value
                                            * distinct values is the number of distinct values divided by the total
                                                number of values (this is inversely related to the quantization error)
                                            * tbp is packet time bandwidth product
                                            * frequency sigma - standard deviation of the frequency
                                            * discontinuity max value - maximum value of discontinuities
                * prediction": (`numpy.array`) - array of a list with a string. Since only concerned with artifacts, these
                                                will all have the following value ['artifact']
            * dropout_fraction (`list`): series of results from the feature metric, starts at first window/ends at last
            * distinct_values_ratio (`list`): series of results from the feature metric, starts at first window/ends at last
            * packet_time_bandwidth_product (`list`): series of results from the feature metric, starts at first window/ends at last
            * frequency_sigma (`list`): series of results from the feature metric, starts at first window/ends at last
            * discontinuity_max_value (`list`): series of results from the feature metric, starts at first window/ends at last
            * artifacts (`list`): series of results from the feature metric, starts at first window/ends at last

    :rtype: dict
    """

    # Call evaluate_stream and read the file in as a stream object
    # Evaluate_stream calls the evaluate function for each trace which evaluates whether it contains signal or an
    # artifact using a random forest ML model
    d = evaluate_stream(
        read(file_path),
        model_path=model_path,
        window_size=window_size,
        stride_length=stride_length,
        database=database_config,
    )

    return d


def evaluate_webservice(
    model_path="/sigpro/qcStatisticsML/models/1621881059_RF-Raw-QC_Balanced.joblib.pkl",
    window_size=60,
    stride_length=50,
    utc_start=None,
    utc_end=None,
    network=None,
    station=None,
    channel=None,
    location=None,
    database_config=None,
):
    """
    Function to evaluate a user-requested SNCL from IRIS webservices, convert it to a stream object, then loop through
    each trace in the stream object (using evaluate_stream function) to evaluate whether it contains signal or an
    artifact using a Random Forest ML model (using evaluate function)

    :param: model_path: training model file path location. Should be relative location as evaluate function will look
                        for full directory path. Default example(
                        '/sigpro/qcStatisticsML/models/1621881059_RF-Raw-QC_Balanced.joblib.pkl')
    :type model_path: str
    :param window_size: rolling window size in seconds
    :type window_size: int
    :param stride_length: length in seconds to increment rolling window
    :type stride_length: int
    :param utc_start: UTC start time for web service request
    :type UTCDateTime object
    :param utc_end: UTC end time for web service request
    :type UTCDateTime object
    :param network: Network code for web service request
    :type: `str`
    :param station: Station code for web service request
    :type: `'str`
    :param channel: Channel code for web service request
    :type: `str`
    :param location: Location code for web service request
    :type: `str`
    :param database: Pycheron database to ingest metrics into
    :type database: pycheron.db.sqllite_db.Database

     returns: d; list of dictionaries containing all the artifacts found for each trace object within
                            respective stream object. The following keys exist:

            * snclq (`str`): station network channel location quality code information for specfic trace
            * start_time (`str`): start time of trace
            * end_time (`str`): end time of trace
            * metric_name (`str`): qcMLMetric string
            * qc_ml_results (`list`): list of dictionaries containing the following key/values
                * window_start: (`UTCDateTime`) - start time of rolling window
                * window_end: (`UTCDateTime`) - end time of rolling window
                * Features: (`numpy array`) - array of feature values: drop out fraction, distinct values ratio, packet time
                                            bandwidth product, frequency sigma, and discontinuity max value.
                                            * dropout fraction is the number of discrete intervals in the trace with N
                                                or more consecutive samples having the same value
                                            * distinct values is the number of distinct values divided by the total
                                                number of values (this is inversely related to the quantization error)
                                            * tbp is packet time bandwidth product
                                            * frequency sigma - standard deviation of the frequency
                                            * discontinuity max value - maximum value of discontinuities
                * prediction": (`numpy.array`) - array of a list with a string. Since only concerned with artifacts, these
                                                will all have the following value ['artifact']
            * dropout_fraction (`list`): series of results from the feature metric, starts at first window/ends at last
            * distinct_values_ratio (`list`): series of results from the feature metric, starts at first window/ends at last
            * packet_time_bandwidth_product (`list`): series of results from the feature metric, starts at first window/ends at last
            * frequency_sigma (`list`): series of results from the feature metric, starts at first window/ends at last
            * discontinuity_max_value (`list`): series of results from the feature metric, starts at first window/ends at last
            * artifacts (`list`): series of results from the feature metric, starts at first window/ends at last

    :rtype: list of dictionaries
    """
    # Set up client
    client = Client("IRIS")
    # Obtain waveforms for given SNCL
    stream = client.get_waveforms(network, station, location, channel, utc_start, utc_end)
    # Loop through traces in stream via evaluate_stream
    # Evaluate_stream calls the evaluate function for each trace which evaluates whether it contains signal or an
    # artifact using a random forest ML model
    d = evaluate_stream(
        stream,
        model_path=model_path,
        window_size=window_size,
        stride_length=stride_length,
        database=database_config,
    )

    return d
