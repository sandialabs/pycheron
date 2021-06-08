**Tutorials**
================

How to run the backend
```````````````````````
There is a config file that must be called in conjunction with `callPycheronMetric.py`. Do the following to run:

`python <path>/<to>/callPycheronMetric.py <path>/<to>/<CONFIG FILE>.yaml`

There is a config file template called `pycheronConfigTemplate.yaml` that contains all of the parameters and their defaults.

For instructions on how to run:

* MSEED files
* Wfdiscs (CSS)
* Obspy Streams
* Wfdisc database (Oracle)


Datasets
-----------------
`callPycheronMetric` accepts five data types:

- Directory (of MSEED files)
- Wfdisc (CSS)
- Obspy Stream
- Singular MSEED file
- Wfdisc database (Oracle)

Depending on the data ingested there are different parameters that can be set. These will be gone over in Section 3.

Directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~


If a directory is being used, pycheron assumes that all files are `.mseed` files. The folder structure should be:

.. code-block:: bash
    
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


With each stream name in the format `<station>_<channel>_<julian day>.mseed`. The data used in this tutorial is already formatted correctly.

Wfdisc (CSS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a Wfidisc (CSS) file is being used the only file that is input into the function is the actually `.wfdisc` file. However, you do need to make sure the `.w` paths listed in the `.wfdisc` file are all relative to the file. For example:

.. code-block:: bash
    
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


and inside the `data.wfdisc`:

.. code-block:: bash

    OWUT   EHZ       1293840000.00000   1927888      774  2011001  1293860963.29000  2096330  100.000000          1.00000         -1.00000 -      o i4 - data_directory               OWUT_EHZ_001.w                            0        -1 10-APR-18        
    OWUT   EHZ       1293861020.31000   1927890      774  2011001  1293861041.29000     2099  100.000000          1.00000         -1.00000 -      o i4 - data_directory               OWUT_EHZ_002.w                      8385320        -1 10-APR-18        
    OWUT   EHZ       1293861041.41000   1927892      774  2011001  1293865831.39000   478999  100.000000          1.00000         -1.00000 -      o i4 - data_directory               OWUT_EHZ_003.w                      8393716        -1 10-APR-18        


Mseed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a single `mseed` file which is read in. It is suggested you use the directory and split the data up if the single stream contains multiple days.

Stream
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is to be used only if inside a python console and not running via terminal. It is assumed that the stream has been either manually loaded using `obspy.read()` or downloaded using the obspy client and then calling the `callPycheronMetric()` as a function


Connecting to a database with Wfdisc tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pycheron currently supports connection to wfdisc databases via SQLAlchemy.


Parameters
------------------------------

The `callPycheronMetric` has over 50 parameters, however the majority of them are all the thresholds for each of the corresponding metrics. What `callPycheronMetric` actually does is combine all of the metrics in the `pycheron.metrics` package and runs them, puts them into a database or csv, and plots figures. It is recommended you throughly read the documentation page on `callPycheronMetric` as it explains every single parameter. 

Running via Command Line
------------------------------

The easiest way to run `callPycheronMetric` is via the command line. Is by using a `config.yaml` file. There is a sample template in the `pycheron` folder that defines all of the default vaules. For this tutorial you can find it in the `data` folder.

Outputs
------------------------------

CSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Output to a CSV file is possible. The output is saved into the following folder structure: 

.. code-block:: bash

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


Database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If desired pycheron can output results into a sqllite database locally on you system. The default name for the database is `pycheron.db` but it can chnaged in the `config.yaml`


callPycheronMetric Example
----------------------------------

The `callPycheronMetric.py` script is installed as a script so you are able to run it directly from the command line

* In the `tutorials/data` directory there is a compressed file that contains some sample data. To de-compress the data move into that directory `cd /tutorial/data` and run 

.. code-block:: bash
    
    tar -xzvf callPycheronMetric_tutorial_data.tar.gz 


* That should create a directory called `pycheron_test` that contains several `mseed` data files. 

* Create a new directory  inside `tutorial/data` called `pycheron_tutorial` by running `mkdir pycheron_tutorial`

* Now open the `callPycheronMetric_tutorial_config.yaml` file in the text editor of your choice. 

* Change the `output_dir` to the absolute path to the `pycheron_tutorial` directory that was recently created. That can be done by navigating to that directory in a terminal and running `pwd` and then copying that output and pasting it as the value for `output_dir`

* You will also need to change the `data` feild. This value will be the absolute path of the directory that has the data that you would like to process. This will be the `pycheron_test` directory. The absolute path can be accessed using the method in the previous step. 

* Change the `datatype` field to be `"dir"` as this will read all `mseed` files in the directory.

* Change `calcAll` to `True` to calculate all of the different metrics or individualy set each metric to either `True` or `False`

You can now run the following command to run the `callPycheronMetric` script

.. code-block:: bash

    python ../pycheron/callPycheronMetric.py data/callPycheronMetric_tutorial_config.yaml


Ouput will look similiar to: 
.. code-block:: bash

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
    



How to run the UI
---------------------------------

Before running the UI, you will need to obtain a [mapbox access token](https://docs.mapbox.com/help/glossary/access-token). After you have created a mapbox account generated a mapbox token, you will need to do the following:

1. Create a file named `.env` in the `pycheron/UI` folder.
2. Create a variable named `MAPBOX_ACCESS_TOKEN` in the `.env` file and assign it the value of your token

An example of what this should look like is provided in `pycheron/UI/.env_example`

After the backend has generated a `.db` file and the corresponding network/station/channel directories with plots, navigate to `pycheron/UI` and run `python createDashUI.py` to generate the UI. Navigate to `localhost:8050` in your browser, and copy/paste the full path to the database in the *CONNECT* box before clicking connect.
