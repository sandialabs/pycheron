## Pycheron Documentation

Pycheron's documentation is created using the `Sphinx AutoDoc` module. To download Sphinx, execute one of the following 
commands in your Python 3.6.9 environment:

### Required Python Packages
The generation of the Sphinx documentation requires, both `sphinx` and `sphinx_bootstrap_theme`. 
If these are not installed in the currnet python environment, they can installed using
* `pip install sphinx`
* `pip install sphinx_bootstap_theme`
* `pip install myst-parser`

### Structure
For every individual Python script (i.e. `correlationMetric.py`) as well as package (i.e. `metrics`) there should be a
corresponding: 

*  `.rst` page located in the `source` folder

In a package `.rst` (i.e. `pycheron.metrics.rst`) page there 
should only be a: 

*  `.. toctree::` containing a list of each individual script`.rst` file, although the `.rst` extension
can be left off


In an individual script (i.e. `pycheron.metrics.correlationMetric.rst`) there should be: 

*  an `.. automodule::` statement followed by the name of the script

AutoDoc will then take documentation from each function's docstring and format it. Thus, when developing a new script or
package, a new `.rst` file will need to be created. 

### Docstring Format

As mentioned above, the Docstring used to document each function is written in `rst` format 
(other formats exist, including `numpy`). A
good cheatsheet can be found [Sphinx tutorial cheatsheet](https://sphinx-tutorial.readthedocs.io/cheatsheet/).

### Figures and Images

Any figures/images that are used in the docstrings or `.rst` pages should be stored in `/source/_static/`.

### LaTeX package
LaTeX is required to build math equations that exist within the documentation. 

At the bottom of the `conf.py` file within the `/docs/source/` directory, the following path declarations are made:

```#----- Options for pngmath ------
# ** Requires LateX Installation **
# Recommend MiKTex for Windows and MacTex for Mac

# For Windows Machines
# pngmath_latex=r"C:\Program Files\MiKTeX 2.9\miktex\bin\x64\latex.exe"
# pngmath_dvipng=r"C:\Program Files\MiKTeX 2.9\miktex\bin\x64\dvipng.exe"

# For Mac
imgmath_latex=r"/Library/TeX/texbin/latex"
imgmath_dvipng=r"/Library/TeX/texbin/dvipng"`
```

Uncomment/comment the correct distribution for your system. The paths included are the standard installation 
locations for the LaTeX distribution. If LaTeX was not installed in the standard installation location, 
file paths will need to be updated. 

### Making the HTML documentation pages

From the `docs` folder, run: 

```
make html
``` 
to create the html documentation pages. Documentation html files will be created and saved into 
`/docs/build/html` directory.

The `index.html` file provides a system diagram of Pycheron as well as: an index to each function within Pycheron; 
a module index, and a search page for finding docs of interest. While the `pycheron.html` provides a table of contents 
of each package and function within Pycheron. 

All html documentation pages provide access to the Pycheron: installation instructions (`Install` tab), 
table of contents of each package and function within Pycheron (`Documentation` tab, same as `pycheron.html` file),
Github repository (`Github` tab), and tutorials (`Tutorials` tab). 