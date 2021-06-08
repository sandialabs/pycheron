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

"""
A script to train a random forest model for use in classifying signals as artifacts or signals (i.e., valid data
without QC issues). To build model run: python train_evaluate_model.py qa_features_labeled.csv <some local directory>

* Code base originally written by Steven Magana-Zook from Lawrence Livermore National Laboratories
  Augmented and adapted for use within Pycheron by Pycheron team
"""

import sys  # Command Line Args
import os  # File path methods
import time  # timestamp filename
import numpy as np
import pandas as pd  # Dataset reading
import sklearn  # Machine learning algorithms
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, accuracy_score, precision_score, recall_score
import imblearn  # Data Augmentation and Balancing
from imblearn import over_sampling


def main(path_to_training_instances_csv, model_output_directory):
    print("Reading in dataset...", end="")
    # Read input file that has features and labels for the model
    ds = pd.read_csv(path_to_training_instances_csv, header=0)
    print("Done")

    # ---- Start code to train a model on balanced data

    # Eliminate the columns that are not features to train on
    # use the top predictive features to make it easier to extract features from new variables
    features = [
        "DROPOUT_FRAC",
        "DISTINCT_VAL_RATIO",
        "TBP",
        "FREQ_SIGMA",
        "DISC_MAX_VALUE",
    ]
    X = ds.loc[:, features]

    # pull out the ground truth labels
    Y = ds.loc[:, "LABEL"]
    label_counts = Y.value_counts()
    print("Before Smote:", label_counts)

    # Applly borderlineSMOTE to balance out the dataset and count the labels
    print("Applying BorderlineSMOTE to balance dataset...", end="")
    borderline_smote = over_sampling.BorderlineSMOTE()
    X, Y = borderline_smote.fit_resample(X, Y)
    print("Done")

    # Turn X, Y back into pandas objects
    X = pd.DataFrame(X, columns=['DROPOUT_FRAC', 'DISTINCT_VAL_RATIO', 'TBP', 'FREQ_SIGMA', 'DISC_MAX_VALUE'])
    Y = pd.Series(Y)

    label_counts = Y.value_counts()
    print("After Smote:\n", label_counts)

    # Split dataset into train-test splits
    print("Shuffling dataset and creating train-test splits...")
    X_train, X_test, y_train, y_test = train_test_split(
        X, Y, test_size=0.25, shuffle=True, random_state=11610
    )

    # Create model based on train splits
    print("Creating RandomForest model based on x_train, y_train...")
    clf_rf = RandomForestClassifier(n_estimators=1000)
    clf_rf.fit(X_train, y_train)
    print("Feature importances: ", list(zip(features, clf_rf.feature_importances_)))

    # Predict test instances
    print("Predicting test instances...")
    x_test_predictions = clf_rf.predict(X_test)

    # Generate metrics on model performance: f1 score, accuracy score, precision score, and recall score
    print("Generating model performance metrics...")
    sklearn_f1_score = f1_score(y_test, x_test_predictions, pos_label="artifact")
    sklearn_acc_score = accuracy_score(y_test, x_test_predictions)
    sklearn_precision_score = precision_score(
        y_test, x_test_predictions, pos_label="artifact"
    )
    sklearn_recall_score = recall_score(
        y_test, x_test_predictions, pos_label="artifact"
    )

    # Print out metric results
    print(("F1: %.2f%%" % (sklearn_f1_score * 100)))
    print(("Accuracy: %.2f%%" % (sklearn_acc_score * 100)))
    print(("Precision: %.2f%%" % (sklearn_precision_score * 100)))
    print(("Recall: %.2f%%" % (sklearn_recall_score * 100)))

    # Create model pickle file. Avoid duplicate filenames and provide natural order to when file created
    model_filename = str(int(time.time())) + "_RF-Raw-QC_Balanced.joblib.pkl"
    model_path = os.path.join(model_output_directory, model_filename)
    print("Saving model to: ", model_path)
    joblib.dump(clf_rf, model_path)


if __name__ == "__main__":
    print("QC Statistic-based Signal/Artifact Classifier")

    print("\tnumpy version: ", np.__version__)
    print("\tpandas version: ", pd.__version__)
    print("\tsklearn version: ", sklearn.__version__)
    print("\timblearn version: ", imblearn.__version__)

    training_instance_path = sys.argv[1]
    model_output_directory = sys.argv[2]

    main(training_instance_path, model_output_directory)
