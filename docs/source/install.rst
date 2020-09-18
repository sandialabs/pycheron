**Install**
-----------

Current Version: **2.1.0**

MasOS and Unix
``````````````
To install:

.. code-block:: bash

   sh build.sh 2.1.0

Windows
```````
To rebuild distribution, run from directory where ``setup.py`` is located:

.. code-block:: bash

   python setup.py sdist bdist_wheel

To install:

.. code-block:: bash

   pip install pycheron-2.1.0-py2-none-any.whl

Known Issues
````````````
- ``obspy 1.1.0`` has been release and updated their module taup, which is used by ``pisces``. An ``ImportError`` my occur. The fix is on line 24 of ``picses util.py``, replace:

   .. code-block:: python

      from obspy.taup import taup

  to

   .. code-block:: python

      from obspy import taup

  A issue has been submited on github to have this bug fixed in ``pisces 0.2.4``. (https://github.com/jkmacc-LANL/pisces/issues/28)


- ``Numpy 1.14.5`` is recommended to use if using Fortran. Higher versions have a bug in ``f2py``

**Running**
-----------

There is a config file that must be called inconjunction with ``callPycheronMetric.py``. When the pycheron lib is installed it adds ``callPycheronMetric.py`` to your ``PATH``,
allowing you to run it directly in the terminal so in order to run in a terminal all that needs to be run is:

.. code-block:: bash

   callPycheronMetric.py <CONFIG FILE>.yaml

There is a config file template called ``pycheronConfigTemplate.yaml`` that contains all of the parameters and their defaults.


Fortran
```````
It is recommended that you use ``gfortran 6.3`` for OSX. Currently, Fortran compilations only work on OSX and Linux, but there is a flag to turn off ``Fortran`` usage.
For Linux and OSX users, the option to use fortran-based scripts is available and recommended for faster runtimes. In order to use fortran a compiler needs to be installed.
It is recommended that `gfortran` be used. For MacOS `gfortran 6.3` is specifically recommended and can be found `here <https://github.com/fxcoudert/gfortran-for-macOS/releases>`_.
For Unix systems gfortran binaries can be found `here <https://gcc.gnu.org/wiki/GFortranBinaries>`_

.. warning:: Fortran option is not available currently in Windows installations. Set the variable ``fortran`` to ``false`` in the config file.

