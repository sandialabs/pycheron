## Pycheron Documentation

Pycheron's documentation is created using `Sphinx AutoDoc` module. Other requirements to make the docs are a LaTeX 
installation to produce the doc's math figures.

### Structure
For every individual python script (i.e. `correlationMetric.py`) as well as package (i.e. `metrics`) there should be a
corresponding `.rst` page located in the `source` folder. In a package `.rst` (i.e. `pycheron.metrics.rst`) page there 
should only be a `.. toctree::` containing a list of each individual script`.rst` files, although the `.rst` extension
can be left off. In a individual script (i.e. `pycheron.metrics.correlationMetric.rst`) there should me an 
`.. automodule::` statement followed by the name of the script. AutoDoc will then take documentation from each
function's docstring and format it. So if you develop a new script or package, you will need to create a new `.rst` file
for those scripts/packages.

### docstring Format

The docstring documenting each function are written in `rst` format (there are other formats, including `numpy`). A
good cheatsheet can be found [here](https://sphinx-tutorial.readthedocs.io/cheatsheet/).

### Figures and Images

Any figures/images that are used in the docstrings or `.rst` pages should be stored in `/source/_static/`

### LaTeX package
LaTeX is required to build the math equations. Inside `conf.py`, at the bottom, the following declarations are made:

```#----- Options for pngmath ------
# ** Requires LateX Installation **
# Recomend MiKTex for Windows and MacTex for Mac

# For Windows Machines
# pngmath_latex=r"C:\Program Files\MiKTeX 2.9\miktex\bin\x64\latex.exe"
# pngmath_dvipng=r"C:\Program Files\MiKTeX 2.9\miktex\bin\x64\dvipng.exe"

# For Mac
imgmath_latex=r"/Library/TeX/texbin/latex"
imgmath_dvipng=r"/Library/TeX/texbin/dvipng"`
```

Uncomment/comment the correct distribution for your system. The paths included are the standard installation spots for
the distribution. Change only if you installed it elsewhere.

### Making the HTML pages

From the `docs` folder, run `make html` to create the html. All of the files will then be created and saved into
 `/docs/build/html`.