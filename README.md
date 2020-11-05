## Pycheron v3.0.0<br>
[![pipeline status](https://gitlab.sandia.gov/lynm/pycheron/badges/tfx-base-python3/pipeline.svg)](https://gitlab.sandia.gov/lynm/pycheron/-/commits/tfx-base-python3) 
[![coverage report](https://gitlab.sandia.gov/lynm/pycheron/badges/tfx-base-python3/coverage.svg)](https://gitlab.sandia.gov/lynm/pycheron/-/commits/tfx-base-python3)


Developed by jbobeck @ SNL and kaaur @ SNL

Python library for seismic data quality control based on IRIS Mustang, IRISSeismic, and seismicRoll

[Dependencies](./setup.py)

To rebuild distribution, run from directory where `setup.py` is located:<br>
`python setup.py sdist bdist_wheel`<br><br>

### Setting up a virtual environment before install
The tool being used for dependency management ([Poetry](https://python-poetry.org/)) will require a virtual environment manager of some sort. A number of options are provided below, though some will work better than others depending on your system architecture (Windows, Mac, or Linux)

* [venv](https://docs.python.org/3/library/venv.html) (comes default with Python 3.3+)
* [PyEnv](https://github.com/pyenv/pyenv) (does not work on Windows)
* [Anaconda](https://docs.anaconda.com/anaconda/install/)
* [PipEnv](https://github.com/pypa/pipenv)

(Alternatively, if you do not use your own environment manager, Poetry will create one for you)

You will need to use a **Python 3.6.9** environment, and **activate it before installing packages and using Pycheron**.

### Installing
Before installing, you will need to download [git-lfs](https://git-lfs.github.com/)

After git-lfs is installed, run the following in the root directory:
```bash
git lfs install
git lfs fetch
git lfs pull
```


The easiest way to install pycheron is to run:

`./build.sh <pycheron-version>`

(current version is `3.0.0`)

This will uninstall any old versions of pycheron and install all the python dependencies.

If you are using Mac OSX, and recieve errors relating to matplotlib, you will 
need to reinstall matplotlib version 2.2.5 and then re-build pycheron using the build script. 

If you encounter `EOFerror` being raised by Pebble, then try installing Pebble version 4.3.10

### Build from scratch

To build the version from scratch run from the top-level directory:

`poetry build`

and then

`cd dist`

from there you can use pip to install the .whl file:

`pip install pycheron-<version>-py3-none-any.whl`

### Know Issues
- `Numpy 1.14.5` is recommended to use if using Fortran. Higher versions have a bug in `f2py`

### Fortran
It is recommended that you use gfortran 6.3 for OSX. Currently Fortran compilations only work on OSX and Linux. If you
 are running windows, there is an option to turn off Fortran in the functions that use it (`staltaMetric`, `psdMetric`,`spikesMetric`, `repeatedAmplitudeMetric`)
by setting the variable `fortran` equal to `False`

## How to run the backend
There is a config file that must be called in conjunction with `callPycheronMetric.py`. When the pycheron lib is installed it adds `callPycheronMetric.py` to your `PATH`,
allowing you to run it directly in the terminal so in order to run in a terminal all that needs to be run is:

`callPycheronMetric.py <CONFIG FILE>.yaml`

There is a config file template called `pycheronConfigTemplate.yaml` that contains all of the parameters and their defaults.

For instructions on how to run:

* MSEED files
* Wfdiscs (CSS)
* Obspy Streams

refer to the [tutorial](./tutorials/callPycheronMetric_tutorial.ipynb).

## How to run the UI
Before running the UI, you will need to obtain a [mapbox access token](https://docs.mapbox.com/help/glossary/access-token). After you have created a mapbox account generated a mapbox token, you will need to do the following:

1. Create a file named `.env` in the `pycheron/UI` folder.
2. Create a variable named `MAPBOX_ACCESS_TOKEN` in the `.env` file and assign it the value of your token

An example of what this should look like is provided in `pycheron/UI/.env_example`

After the backend has generated a `.db` file and the corresponding network/station/channel directories with plots, navigate to `pycheron/UI` and run `python createDashUI.py` to generate the UI. Navigate to `localhost:8050` in your browser, and copy/paste the full path to the database in the *CONNECT* box before clicking connect.

## Formatting
For any changes, code must follow [Black](https://github.com/psf/black) formatting rules. You can install/run Black against your changes before pushing with:
* `pip install black`
* `black <new-or-changed-py-files>`

## Documentation
Refer to ["Making the HTML pages" in the docs folder](./docs/README.md#making-the-html-pages) to create Sphinx documentation.

