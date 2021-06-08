## Pycheron v3.0.0<br>
[![pipeline status](https://gitlab.sandia.gov/lynm/pycheron/badges/tfx-base-python3/pipeline.svg)](https://gitlab.sandia.gov/lynm/pycheron/-/commits/tfx-base-python3) 
[![coverage report](https://gitlab.sandia.gov/lynm/pycheron/badges/tfx-base-python3/coverage.svg)](https://gitlab.sandia.gov/lynm/pycheron/-/commits/tfx-base-python3)


Developed by Kale Aur (kaaur@sandia.gov), Jessica Bobeck (jbobeck@sandia.gov), Anthony Alberti (aalber@sandia.gov), 
and Phillip Kay (prkay@sandia.gov)

Python library for seismic data quality control based on IRIS's IRISMustangMetrics, IRISSeismic, and seismicRoll R 
packages:

*   IRISMustangMetrics R Cran Package
    (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, March 20). CRAN-Package IRISMustangMetrics.
    The Comprehensive R Archive Network. Retrieved from
    https://cran.r-project.org/web/packages/IRISMustangMetrics/index.html)

*   IRISSeismic R Cran Package
    (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
    The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index)
    
*   Code originally ported from seismicRoll R Cran Package
    (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive R
    Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index)

*   IRISMustangMetrics R Cran Package
    (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, March 20). CRAN-Package IRISMustangMetrics.
    The Comprehensive R Archive Network. Retrieved from
    https://cran.r-project.org/web/packages/IRISMustangMetrics/index.html)

*   IRISSeismic R Cran Package
    (Callahan, J., R. Casey, G. Sharer, M. Templeton, and C. Trabant (2019, Oct 22). CRAN-Package IRISSeismic.
    The Comprehensive R Archive Network. Retrieved from https://cran.r-project.org/web/packages/IRISSeismic.index)
    
*   Code originally ported from seismicRoll R Cran Package
    (Callahan, J., R. Casey, M. Templeton, and G. Sharer (2020, July 8). CRAN-Package seismicRoll. The Comprehensive R
    Archive Network. Retrieved from https://cran.r-project.org/web/packages/seismicRoll.index)

*   Many of the baseline function thresholds are based on IRIS's thresholds found in the section labeled 
    `Some Metric Tests for Broadband Data`: 
    IRIS: Tutorials: Seismic Data Quality Assurance Using IRIS MUSTANG Metrics. (2016, April 26). Tutorials: Seismic Data Quality Assurance Using IRIS MUSTANG 
    Metrics. https://ds.iris.edu/ds/nodes/dmc/tutorials/seismic-data-quality-assurance-using-iris-mustang-metrics/


### Setting up a virtual environment before install
The tool used for dependency management ([Poetry](https://python-poetry.org/)) requires the use of a virtual environment manager. 
Below are a few virtual environment management options; however, depending on your system architecture (i.e., Windows, Mac, or Linux), some may work better than others:

* [venv](https://docs.python.org/3/library/venv.html) (comes default with Python 3.3+)
* [PyEnv](https://github.com/pyenv/pyenv) (does not work on Windows)
* [Anaconda](https://docs.anaconda.com/anaconda/install/)
* [PipEnv](https://github.com/pypa/pipenv)

Alternatively, if you do not use your own environment manager, Poetry will create one for you. 

Note, a **Python 3.6.9** environment is required and must be *activated before installing packages and using Pycheron. 

### Installing
Before installing, you will need to download [poetry](https://python-poetry.org) and [git-lfs](https://git-lfs.github.com/)

To download poetry, execute one of the following commands in your Python 3.6.9 environment:

```
pip install poetry
```

Or if using Anaconda:

```
conda install poetry 
```

After git-lfs is installed, run the following commands in the Pycheron root directory:
```bash
git lfs install
git lfs fetch
git lfs pull
```


The easiest way to install Pycheron is to run:

`./build.sh 3.0.0`

where 3.0.0 is the version number of Pycheron. This bash script will uninstall any old versions of Pycheron and install all necessary Python dependencies.

If you encounter an `EOFerror` raised by Pebble, then try installing Pebble version 4.3.10 and re-run the `build.sh` script.

### Build from scratch

To build the version from scratch run from the Pycheron top-level directory:

`poetry build`

and then:

`cd dist`

Then use pip to install the .whl file:

`pip install pycheron-<version>-py3-none-any.whl`

### Installing Oracle dependencies
`cx_Oracle` requires that Oracle Instant Client be installed on your system. Instructions for installation can be found [here](https://cx-oracle.readthedocs.io/en/latest/user_guide/installation.html#install-oracle-client). Be sure to export the `LD_LIBRARY_PATH` environment variable before using.

### Known Issues
If you are using Mac OSX with an Anaconda environment and receive errors relating to matplotlib, you will 
need to reinstall matplotlib version 2.2.5 using `conda install matplotlib=2.2.5" and re-build Pycheron using the build script. 

### Fortran
It is recommended that you use gfortran 6.3 for OSX. Currently Fortran compilations only work on OSX and Linux. If you
are using a Windows operating system, there is an option to turn off Fortran in the functions that use it (e.g., `staltaMetric`, `psdMetric`,`spikesMetric`, `repeatedAmplitudeMetric`)
by setting the variable `fortran` equal to `False`

## How to run Pycheron's backend
Data is ingested into Pycheron by directing the package to a local directory 
containing data files, a single data file, a CSS3.0 wfdisc table (flatfile), an ObsPy stream, or a database containing 
CSS3.0 formatted tables. Pycheron’s main interface is via a wrapper code, `callPycheronMetric.py`, that reads from a 
configuration YAML Ain’t Markup Language (YAML) file, `pycheronConfigTemplate.yaml`.   
<br> 
This template file specifies: 

   1) from where Pycheron reads input data; 
   2) what QC metrics are calculated; 
   3) QC metric parameter and threshold settings; 
   4) plotting information; and 
   5) various other parametrizations.

`callPycheronMetric` combines all of the metrics in the `pycheron.metrics` package and executes them. For each QC issue, 
weighting is assigned based on user input and an overall QC summary is created on a per station basis. 
All results are stored in a Sqlite3 database and can either be visualized in a user interface (UI) via Plotly’s 
Dashboard Web UI, Dash, or output to comma-separated value (CSV) files that are more readily processed through an 
automated pipeline. 

`pycheronConfigTemplate.yaml` is a configuration template containing all configurable input parameters 
and their default values. The `pycheronConfigTemplate.yaml` file can be found in the top-level directory of Pycheron.

Please refer to the sections below for instructions on how to run Pycheron with the following data types:

* Directory (of MSEED files)
* Wfdiscs (CSS)
* Obspy Streams
* Singular MSEED file 
* Wfdisc database (Oracle)

## Datasets

Pycheron can ingest several different data formats. `callPycheronMetric` currently accepts five data types:

- Directory (of MSEED files) (datatype = `dir`)
- Wfdisc (CSS) (datatype = `wfdisc`)
- Obspy Stream (datatype = `stream`)
- Singular MSEED file (datatype = `mseed`)
- Wfdisc database (Oracle) (datatype = `wfdb`)

All data formats are converted into ObsPy Stream or Trace objects for use inside Pycheron. 
Thus, any format read by ObsPy can also be used within Pycheron. The data type needs to be defined in the YAML
configuration file, as does the data file's directory path. 
Depending on the data ingested there are different parameters that can be set; these will be discussed in the 
`callPycheronMetric Example` section. 

### Directory (datatype = `dir`)

If a directory is utilized, Pycheron assumes that all files are `MSEED` formatted files. The following 
folder structure is required:

```
My Folder
    BGU_EHE_001.mseed
    BGU_EHN_001.mseed
    BGU_EHE_002.mseed
    BGU_EHN_002.mseed
    JPU_EHE_001.mseed
    JPU_EHN_001.mseed
    JPU_EHE_002.mseed
    JPU_EHN_002.mseed
    TCRU_EHE_001.mseed
    TCRU_EHN_001.mseed
    TCRU_EHE_002.mseed
    TCRU_EHN_002.mseed
```

Each file should be formatted as follows: `<station>_<channel>_<julian day>.mseed`. 
The data utilized within this tutorial is already formatted correctly. 
If the files within the directory are not formatted correctly, 

### Wfdisc (CSS) (datatype = `wfdisc`)
If a Wfdisc (CSS) flat file is utilized, the `.wfdisc` file is the only input file needed.  
However, it is important to ensure that the `.w` file paths in the dir column listed in the `.wfdisc` file are 
relative to the file. 

For example, using the following directory structure below:

```
My Folder
    data.wfdisc
    data_directory
        OWUT_EHZ_001.w 
        OWUT_EHZ_002.w 
        OWUT_EHZ_003.w 
        .
        .
        .
        waveformN.w
```

The dir column in the `data.wfdisc` file should appropriately point to the the My Folder/data_directory for each `.w` file 
so that Pycheron is able to properly retrieve the respective file:

```
OWUT   EHZ       1293840000.00000   1927888      774  2011001  1293860963.29000  2096330  100.000000          1.00000         -1.00000 -      o i4 - /My Folder/data_directory/               OWUT_EHZ_001.w                            0        -1 10-APR-18        
OWUT   EHZ       1293861020.31000   1927890      774  2011001  1293861041.29000     2099  100.000000          1.00000         -1.00000 -      o i4 - /My Folder/data_directory/               OWUT_EHZ_002.w                      8385320        -1 10-APR-18        
OWUT   EHZ       1293861041.41000   1927892      774  2011001  1293865831.39000   478999  100.000000          1.00000         -1.00000 -      o i4 - /My Folder/data_directory/               OWUT_EHZ_003.w                      8393716        -1 10-APR-18        
```

### Mseed (datatype = `mseed`)

This datatype indicates that the user wishes to read in a single `mseed` formatted file.  
If the file contains multiple days worth of data, it is recommended to use the directory datatype and split
the data into smaller segments. 

### Stream (datatype = `stream`)

This datatype should only be used when inside a Python console and not while executing via the command line. 
When this type is specified, it is assumed that the Stream has either been  manually loaded using 
`obspy.read()` or downloaded using the ObsPy client and then calling `callPycheronMetric()` as a function.

### Connecting to a database with Wfdisc tables (datatype = `wfdb`)

Pycheron currently supports connection to Oracle databases via SQLAlchemy. This datatype option supports 
reading Wfdisc database tables. 

### `callPycheronMetric` Input Parameters

The `callPycheronMetric.py` script is installed as a script so that it can be executed directly from the command line. 
The `callPycheronMetric` has over 50 input parameters specified, with a majority of them related to setting default 
threshold values for each of the corresponding metrics. These thresholds determine whether a QC issue is flagged. 
Thus, it is recommended to thoroughly read through the documentation page within the  `callPycheronMetric` script as 
each parameter is well-documented. 

### Running Pycheron via Command Line

The easiest way to run `callPycheronMetric` is via the command line with the `pycheronConfigTemplate.yaml` file. 
After installation and configuration are complete, Pycheron can be run via the terminal 
with the following command:

`python <path>/<to>/callPycheronMetric.py <path>/<to>/<CONFIG FILE>.yaml`

For this tutorial, the `pycheronConfigTemplate.yaml` is accessible via the `data` folder.

## Pycheron Output Options

### CSV output 

Output to a CSV file is possible. The output is saved into the following folder structure: 

```
<output_directory>
    <network_dir>
        <station_dir>
            <channel_1_dir>
                <metric>.csv
                <metric_plot>.png
            <channel_2_dir>
                <metric>.csv
                <metric_plot>.png
            <station_level_metrics>.csv
            <station_level_plot>.png
         <network_level_metric>.csv
         <network_level_plot>.png
```

### Database output 

If desired, Pycheron can output results into a local sqllite database. 
The default name for the database is `pycheron.db` but this can be changed within the 
`pycheronConfigTemplate.yaml` file.  


##  callPycheronMetric Examples

This section will walk users through how to use Pycheron with each input data type. 

Sample data files for this tutorial live withIn the `tutorials/data` directory. 

To process the data, execute the following steps: 

1. To de-compress the data, first cd into the aforementioned directory: 

   `cd /tutorials/data` 

    then run:

    ```
    tar -xzvf callPycheronMetric_tutorial_data.tar.gz 
    ```

    This should create a directory called `pycheron_test` that contains several `mseed` data files, spanning 
    2011/01/01-2011/01/06, for 5 stations within the University of Utah (UU) seismic network: BGU, CTU, HVU, MTPU, and 
    ZNPU. 

2. Create a new directory inside `tutorials/data` called `pycheron_tutorial` by running the following command:

   `mkdir pycheron_tutorial`

3. Open the `callPycheronMetric_tutorial_config.yaml` file in a text editor of your choice. 

   3.1. Change the `output_dir` in the `callPycheronMetric_tutorial_config.yaml` to the absolute path of the 
        `pycheron_tutorial` directory that was recently created above. One way to obtain the absolute path is
        to navigate to `pycheron_tutorial` directory in a terminal and execute the `pwd` command and then copy 
        and paste it as the value for the `output_dir` input parameter:
        
        `output_dir = "/Users/username/pycheron/tutorials/data/pycheron_tutorial"`

   3.2. Next, change the `data` field to the absolute path of the directory that has the data that you would 
        like to process. For this tutorial, the `data` field should point to the `pycheron_test` directory:
        
        `data = "/Users/username/pycheron/tutorials/data/pycheron_test"`

   3.3. Change the `datatype` field to be `"dir"`. This setting will read all `mseed` files within the directory:
   
        `datatype = "dir"`

   3.4. Change `calcAll` to `True` to calculate all of the different metrics or individually set each metric to 
        either `True` or `False` to delineate which metrics to calculate:
        
        `calcAll = True`
        
   3.5. For simplicity within this tutorial, it is recommended to keep all other default thresholds and input parameters.
        However, if the user would like to experiment with updating other input parameters it is recommended to
        thoroughly read through the documentation to learn which input parameters map to which metric and what their 
        default settings are. Remember, if on a Windows system, to turn off Fortran in the functions that use it. 

4. Execute the following command to run the `callPycheronMetric` script: 

    ```
    python ../../pycheron/callPycheronMetric.py callPycheronMetric_tutorial_config.yaml
    ```

Example output while Pycheron is processing data: 
```
-----------------------------------------------
Plotting UU.BGU Jul Date: 002
-----------------------------------------------
/Users/prkay/Workspace/gh-pycheron/tutorials/data/pycheron_tutorial/pycheron.db already exists. Connecting to /Users/prkay/Workspace/gh-pycheron/tutorials/data/pycheron_tutorial/pycheron.db...
No results found
No results found
Finished PDFgrid and Line plots: UU.BGU
/Users/prkay/Workspace/gh-pycheron/tutorials/data/pycheron_tutorial/pycheron.db already exists. Connecting to /Users/prkay/Workspace/gh-pycheron/tutorials/data/pycheron_tutorial/pycheron.db...
Finished stationNoisePlot: UU.BGU
/Users/prkay/Workspace/gh-pycheron/tutorials/data/pycheron_tutorial/pycheron.db already exists. Connecting to /Users/prkay/Workspace/gh-pycheron/tutorials/data/pycheron_tutorial/pycheron.db...
Finished psdPlot: UU.BGU
/Users/prkay/Workspace/gh-pycheron/tutorials/data/pycheron_tutorial/pycheron.db already exists. Connecting to /Users/prkay/Workspace/gh-pycheron/tutorials/data/pycheron_tutorial/pycheron.db...
Finished pdfPlot: UU.BGU
....
....

```

## How to run the Pycheron User Interface (UI)
Before running the UI, users will need to obtain a [mapbox access token](https://docs.mapbox.com/help/glossary/access-token). 
After a mapbox account and a mapbox token have been created, execute the following steps:

1. Create a file named `.env` in the `pycheron/UI` folder.
2. Create a variable named `MAPBOX_ACCESS_TOKEN` in the `.env` file and assign it the value of the mapbox token created 
   above

An example of what this should look like is provided in `pycheron/UI/.env_example` file.

After a database has been created (i.e., a `.db` file exists) and the corresponding network/station/channel directories 
with plots exist: 

1.  Navigate to `pycheron/UI` folder and execute the following command in the terminal: 

    `python createDashUI.py` to generate the UI. 

2.  Navigate to `localhost:8050` in your browser

3.  Copy/paste the full path to the database in the *CONNECT* box before clicking `connect`.

## Formatting
For any changes, code must follow [Black](https://github.com/psf/black) formatting rules. Users can install/run Black 
against any changes before pushing with:
* `pip install black`
* `black <new-or-changed-py-files>`

## Developing Locally
For users who intend to make changes to Pycheron before running, it is suggested to execute the following commands
so that Pycheron will not have to be rebuilt with each change:  

* From the root directory, run `cd dist` after Pycheron has been built
* Extract the tarball with `tar -xzvf pycheron-<verion-number-here>.tar.gz`
* Run `cd pycheron-<version-number-here>`, and copy/paste the `setup.py` into the Pycheron root directory 
* From the root directory, run `pip install -e .`

## Documentation
Refer to ["Making the HTML pages" in the docs folder](./docs/README.md#making-the-html-pages) tutorial to create 
Pycheron's Sphinx documentation. 

