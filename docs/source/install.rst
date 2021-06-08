**Install**
-----------

Current Version: **3.0.0**

Setting up a virtual environment before install
````````````````````````````````````````````````````````
The tool being used for dependency management ([Poetry](https://python-poetry.org/)) will require a virtual environment manager of some sort. A number of options are provided below, though some will work better than others depending on your system architecture (Windows, Mac, or Linux)

* [venv](https://docs.python.org/3/library/venv.html) (comes default with Python 3.3+)
* [PyEnv](https://github.com/pyenv/pyenv) (does not work on Windows)
* [Anaconda](https://docs.anaconda.com/anaconda/install/)
* [PipEnv](https://github.com/pypa/pipenv)

(Alternatively, if you do not use your own environment manager, Poetry will create one for you)

You will need to use a **Python 3.6.9** environment, and **activate it before installing packages and using Pycheron**.


Installing
``````````````````````
Before installing, you will need to download [git-lfs](https://git-lfs.github.com/)

After git-lfs is installed, run the following in the root directory:
```bash
git lfs install
git lfs fetch
git lfs pull
```

The easiest way to install pycheron is to run:

`./build.sh`

This will uninstall any old versions of pycheron and install all the python dependencies.

If you encounter `EOFerror` being raised by Pebble, then try installing Pebble version 4.3.10

Build From Scratch
`````````````````````````````

To build the version from scratch run from the top-level directory:

`poetry build`

and then

`cd dist`

from there you can use pip to install the .whl file:

`pip install pycheron-<version>-py3-none-any.whl`

Installing Oracle dependencies
``````````````````````````````````````````

`cx_Oracle` requires that Oracle Instant Client be installed on your system. Instructions for installation can be found [here](https://cx-oracle.readthedocs.io/en/latest/user_guide/installation.html#install-oracle-client). Be sure to export the `LD_LIBRARY_PATH` environment variable before using.


Known Issues
``````````````````
If you are using Mac OSX with an Anaconda environment and recieve errors relating to matplotlib, you will 
need to reinstall matplotlib version 2.2.5 using `conda install matplotlib=2.2.5" and re-build pycheron using the build script. 


- ``Numpy 1.14.5`` is recommended to use if using Fortran. Higher versions have a bug in ``f2py``

**Running**
-----------

There is a config file that must be called inconjunction with ``callPycheronMetric.py``. When the pycheron lib is installed it adds ``callPycheronMetric.py`` to your ``PATH``,
allowing you to run it directly in the terminal so in order to run in a terminal all that needs to be run is:

.. code-block:: bash

   callPycheronMetric.py <CONFIG FILE>.yaml

There is a config file template called ``pycheronConfigTemplate.yaml`` that contains all of the parameters and their defaults.


Fortran
`````````````````
It is recommended that you use ``gfortran 6.3`` for OSX. Currently, Fortran compilations only work on OSX and Linux, but there is a flag to turn off ``Fortran`` usage.
For Linux and OSX users, the option to use fortran-based scripts is available and recommended for faster runtimes. In order to use fortran a compiler needs to be installed.
It is recommended that `gfortran` be used. For MacOS `gfortran 6.3` is specifically recommended and can be found `here <https://github.com/fxcoudert/gfortran-for-macOS/releases>`_.
For Unix systems gfortran binaries can be found `here <https://gcc.gnu.org/wiki/GFortranBinaries>`_

.. warning:: Fortran option is not available currently in Windows installations. Set the variable ``fortran`` to ``false`` in the config file.


