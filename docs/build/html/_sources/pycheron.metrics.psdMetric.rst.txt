pycheron.metrics.psdMetric
--------------------------

Performs a spectral analysis on seismic traces within stream object and returns PSD metrics. This metric returns a
list of average power spectra associated with a stream's data by breaking the stream into
chunks and calculating the spectrum for each. Methodology based on the 2004 McNamara and Buland paper. [#f1]_

.. automodule:: pycheron.metrics.psdMetric
    :members:
    :show-inheritance:

.. rubric:: Footnotes

.. [#f1] McNamara and Buland. 2004. "Ambient Noise Levels in the Continental United States". https://profile.usgs.gov/myscience/upload_folder/ci2012Feb2217152844121McNamaraBuland_BSSA.pdf