"""Simple workflow for using the Bioconductor package `markeR` from Python via rpy2.

source .venv/bin/activate
pip install -r requirements.txt
python markeR_to_python.py --tutorial   # now should execute successfully

This script shows how to:

1. Configure rpy2 and ensure an R environment is available.
2. Install Bioconductor and markeR if not already installed.
3. Load R functions into Python.
4. Demonstrate a markeR analysis using example data.

Notes:
- You must have R (>=4.5) installed on your system.
- Install the Python package `rpy2` in the same environment where this script runs:
    pip install rpy2

For more details on markeR see https://bioconductor.org/packages/markeR

This module also includes utilities (`plot_r_expression`, `plot_r_function`)
that open an R graphics device and capture plots as PNG files.  When run inside a
Jupyter notebook the plots are automatically displayed inline; otherwise the
images are saved to a temporary file whose path is printed.
"""

from __future__ import annotations

import sys

# check that required Python libraries are installed before proceeding
_missing = []
try:
    import numpy  # used by rpy2 and examples
except ImportError:  # pragma: no cover - dependency check
    _missing.append("numpy")
try:
    import pandas  # examples use it for conversions
except ImportError:  # pragma: no cover
    _missing.append("pandas")
try:
    import rpy2  # primary bridge to R
except ImportError:  # pragma: no cover
    _missing.append("rpy2")

if _missing:
    sys.exit(
        "The following Python packages are required but not installed: %s.\n"
        "Please install them (e.g. `pip install -r requirements.txt`)." %
        ", ".join(_missing)
    )

# rpy2 imports
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects import conversion
from rpy2.robjects.packages import importr, isinstalled

# note: pandas2ri.activate() is deprecated; we use conversion contexts when
# converting.  helpers below wrap the recommended API.


def _to_py(obj):
    """Convert an R object to a pandas/numpy equivalent."""
    with conversion.localconverter(ro.default_converter + pandas2ri.converter):
        return conversion.rpy2py(obj)

# utilities for inline plotting (e.g. in Jupyter notebooks)
import os
import tempfile
try:
    from IPython.display import Image, display
    _HAS_IPYTHON = True
except ImportError:  # not running in notebook
    _HAS_IPYTHON = False


def _r_open_png(width=800, height=600, filename=None):
    """Start an R PNG device, returning the filename used."""
    if filename is None:
        filename = tempfile.mktemp(suffix=".png")
    else:
        # Convert to absolute path to ensure R saves to the intended location
        filename = os.path.abspath(filename)
    # ensure the path exists
    dirname = os.path.dirname(filename)
    if dirname and not os.path.isdir(dirname):
        os.makedirs(dirname, exist_ok=True)
    ro.r(f'png("{filename}", width={width}, height={height})')
    return filename


def _r_close_device():
    """Close the active R graphics device."""
    ro.r('dev.off()')


def plot_r_expression(expr: str, width=800, height=600, filename=None, display_plot=True):
    """Evaluate an R expression that produces a plot and optionally display it.

    Parameters
    ----------
    expr : str
        R code string that generates a plot when evaluated.
    width, height : int
        Dimensions for the PNG device in pixels.
    filename : str or None
        Path to which the image should be saved. If None a temporary file
        will be created.
    display_plot : bool
        If True and running under IPython, display the resulting PNG
        inline. Otherwise the path is printed.

    Returns
    -------
    str
        The path to the saved PNG file.
    """
    fname = _r_open_png(width=width, height=height, filename=filename)
    ro.r(expr)
    _r_close_device()
    if display_plot and _HAS_IPYTHON:
        display(Image(filename=fname))
    else:
        print(f"plot written to {fname}")
    return fname


def plot_r_function(func_name: str, *args, width=800, height=600, filename=None,
                    display_plot=True, **kwargs):
    """Call an R plotting function by name and save/display result.

    Any positional and keyword arguments are converted to their R
    equivalents by rpy2.

    Example::

        plot_r_function('PlotScores', data=counts_example, metadata=metadata_example,
                        gene_sets=genesets_example)
    """
    fname = _r_open_png(width=width, height=height, filename=filename)
    rfunc = ro.r[func_name]
    # convert kwargs into ro objects (rpy2 handles this automatically)
    rfunc(*args, **kwargs)
    _r_close_device()
    
    # Try to display if requested; if in actual IPython/Jupyter context, display inline;
    # otherwise just print the path.
    if display_plot and _HAS_IPYTHON:
        try:
            # Check if we're actually in an interactive IPython shell (not just that it's installed)
            from IPython import get_ipython
            ipython = get_ipython()
            if ipython is not None:
                display(Image(filename=fname))
            else:
                print(f"plot written to {fname}")
        except Exception:
            # If anything goes wrong (file not ready, not in IPython context, etc.)
            print(f"plot written to {fname}")
    else:
        print(f"plot written to {fname}")
    return fname


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def ensure_bioc_installed() -> None:
    """Install Bioconductor's package manager if it is not already present."""
    biocinstaller = "BiocManager"
    if not isinstalled(biocinstaller):
        ro.r('install.packages("{0}")'.format(biocinstaller))
    ro.r('suppressMessages(require({0}))'.format(biocinstaller))


def install_markeR() -> None:
    """Install the markeR package from Bioconductor if not already installed.

    Uses BiocManager to perform the installation.  After running this function
    the package should be loadable via `importr("markeR")`.
    """
    ensure_bioc_installed()
    if not isinstalled("markeR"):
        ro.r('BiocManager::install("markeR", ask=FALSE, update=FALSE)')
    ro.r('library(markeR)')


def get_markeR_functions() -> ro.Environment:
    """Return the markeR namespace so that functions can be accessed conveniently.

    Example:
        mark = get_markeR_functions()
        scores = mark.CalculateScores(data=counts, metadata=metadata, gene_sets=genesets, method="logmedian")
    """
    install_markeR()
    # importing via importr is more reliable than accessing `ro.r['markeR']`.
    try:
        return importr("markeR")
    except Exception as e:
        raise RuntimeError("Unable to load markeR package: %s" % e)





# ---------------------------------------------------------------------------
# Tutorial helpers using markeR example data
# ---------------------------------------------------------------------------

def load_benchmark_examples():
    """Load the built-in example data and gene sets from the markeR package.

    Returns a tuple `(counts, metadata, genesets)` where each element is an R
    object.  You can convert them to pandas objects if desired.
    """
    # ensure package is installed and loaded
    install_markeR()
    # load the three example datasets provided by the vignette
    ro.r('data("genesets_example", package="markeR")')
    ro.r('data("counts_example", package="markeR")')
    ro.r('data("metadata_example", package="markeR")')
    genesets = ro.r('genesets_example')
    counts = ro.r('counts_example')
    metadata = ro.r('metadata_example')
    return counts, metadata, genesets


def tutorial_benchmark(output_file=None):
    """Demonstrate a small benchmarking mode example from the markeR vignette.
    
    Parameters
    ----------
    output_file : str or None
        If provided, saves the display output to a file using the R graphics device.
        Note: The markeR::PlotScores function outputs to the active graphics device.
    """
    print("-- loading example data from markeR")
    counts, metadata, genesets = load_benchmark_examples()

    # show dimensions of the data
    print("counts matrix dimensions:", ro.r('dim')(counts))
    print("metadata dimensions:", ro.r('dim')(metadata))
    print("available gene sets:", list(genesets.names))

    # run CalculateScores (logmedian method) as in the tutorial
    calculate = ro.r['CalculateScores']
    print("-- calculating scores using logmedian method")
    df_scores = calculate(data=counts,
                          metadata=metadata,
                          method="logmedian",
                          gene_sets=genesets)

    # df_scores is an R list with one element per gene set; convert first one
    # to pandas for display
    first_name = list(df_scores.names)[0]
    r_first = df_scores.rx2(first_name)
    try:
        import pandas as pd
        pd_first = _to_py(r_first)
        print(f"first gene set ({first_name}) scores (head):\n", pd_first.head())
    except ImportError:
        print("pandas not available; skipping conversion of results to DataFrame")

    # Generate plot using PlotScores
    # Note: markeR's PlotScores function creates an interactive plot or writes to the current device
    print("-- generating a simple score plot")
    
    if output_file:
        # Set up PNG device
        output_file = os.path.abspath(output_file)
        dirname = os.path.dirname(output_file)
        if dirname and not os.path.isdir(dirname):
            os.makedirs(dirname, exist_ok=True)
        ro.r(f'png("{output_file}", width=800, height=400)')
        print(f"   saving to: {output_file}")
    
    # Call the plotting function and force evaluation of the returned plot
    ro.r('''
    p <- PlotScores(
        data = counts_example,
        metadata = metadata_example,
        gene_sets = genesets_example,
        Variable = "Condition",
        method = "logmedian",
        nrow=1
    )
    print(p)
    ''')
    
    if output_file:
        # Close device
        ro.r('dev.off()')
        print(f"plot saved to {output_file}")

# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    output_file = None
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        if sys.argv[1] == "--tutorial":
            # Check if --output flag is present
            if len(sys.argv) > 2 and sys.argv[2].startswith("--output"):
                if sys.argv[2] == "--output" and len(sys.argv) > 3:
                    output_file = sys.argv[3]
                elif "=" in sys.argv[2]:
                    output_file = sys.argv[2].split("=", 1)[1]
            tutorial_benchmark(output_file=output_file)
        else:
            print("usage: python markeR_to_python.py --tutorial [--output FILENAME]")
            print("  --tutorial              : load markeR example data and compute logmedian scores")
            print("  --output FILENAME       : save plot to specified PNG file (optional)")
            print("                            example: python markeR_to_python.py --tutorial --output my_plot.png")
    else:
        print("usage: python markeR_to_python.py --tutorial [--output FILENAME]")
        print("See the module docstring for more details.")