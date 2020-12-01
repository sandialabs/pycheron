# -*- coding: utf-8 -*-

# DO NOT EDIT THIS FILE!
# This file has been autogenerated by dephell <3
# https://github.com/dephell/dephell

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

readme = ""

setup(
    long_description=readme,
    name="pycheron",
    version="3.0.0",
    python_requires="==3.6.9",
    author="Your Name",
    author_email="you@example.com",
    packages=[
        "pycheron",
        "pycheron.UI",
        "pycheron.dataAcq",
        "pycheron.db",
        "pycheron.metricStore",
        "pycheron.metrics",
        "pycheron.metrics.tests",
        "pycheron.plotting",
        "pycheron.psd",
        "pycheron.psd.noise",
        "pycheron.rollseis",
        "pycheron.sigpro",
        "pycheron.test",
        "pycheron.test_callPycheronMetric",
        "pycheron.util",
    ],
    package_dir={"": "."},
    package_data={
        "pycheron": [
            "test_data/*.HHZ",
            "test_data/*.MSEED",
            "test_data/*.csv",
            "test_data/*.db",
            "test_data/*.mseed",
            "test_data/*.site",
            "test_data/*.sitechan",
            "test_data/*.sql",
            "test_data/*.w",
            "test_data/*.wfdisc",
            "test_data/test_data/*.mseed",
        ],
        "pycheron.UI": [
            "*.png",
            "additional_default_css_templates_dash/*.css",
            "assets/*.css",
            "css/*.css",
            "images/*.png",
            "js/*.js",
            "js/hold/*.js",
        ],
        "pycheron.metricStore": ["data/000/*.h5"],
        "pycheron.metrics": [
            "*.mseed",
            "calibration_ml_objects/*.joblib",
            "calibration_ml_objects/*.pkl",
            "repAmps/*.dat",
            "repAmps/*.f90",
        ],
        "pycheron.plotting": ["*.txt"],
        "pycheron.psd": ["fftpack/*.dat", "fftpack/*.f", "fftpack/*.pyf"],
        "pycheron.psd.noise": ["outliers/*.f"],
        "pycheron.rollseis": ["*.f90"],
        "pycheron.test": ["test_data/*.mseed", "test_data/*.yml"],
        "pycheron.test_callPycheronMetric": ["*.yaml"],
    },
    install_requires=[
        "certifi==2019.3.9",
        "cmdline-ispaq",
        "coverage==4.5.4",
        "cx-oracle==7.3.0",
        "dash==0.43.0",
        "joblib==0.16.0",
        "llvmlite==0.33.0",
        "matplotlib==2.2.5",
        "mock==3.0.5",
        "numba==0.51.1",
        "numpy==1.19.1",
        "obspy==1.2.2",
        "pandas==0.24.2",
        "pathlib2==2.3.3",
        "pebble==4.3.8",
        "pisces==0.3.2",
        "pyhht==0.1.0",
        "pytest==4.6.2",
        "pytest-cov==2.8.1",
        "python-dotenv==0.14.0",
        "pywavelets==1.0.2",
        "pyyaml==5.1",
        "scikit-learn==0.20.3",
        "scipy==1.2.1",
        "snakeviz==2.1.0",
        "sqlalchemy==1.3.17",
        "statsmodels==0.9.0",
        "tables==3.5.1",
        "testfixtures==6.14.1",
    ],
    dependency_links=[
        "git+https://github.com/iris-edu/ispaq.git@master#egg=cmdline-ispaq"
    ],
)
