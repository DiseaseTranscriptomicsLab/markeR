# markeR Python Bridge

This workspace provides two simple Python helpers for using the
Bioconductor R package **markeR** via `rpy2`.

- `markeR_to_python.py` provides a minimal example workflow 
  demonstrating how to load the example data and call a specific 
  function, serving as a template for users who wish to structure 
  their own analysis script in Python.
- `run_marker_function.py` is a generic command‑line wrapper that can
  invoke *any* function exported by the markeR package.

## Prerequisites

* R (>=4.5) installed and on your `PATH`.
* A Python virtual environment.  A `requirements.txt` file is provided in
  this folder listing the needed packages (`rpy2`, `pandas`, `numpy` plus
  optional `ipython`/`jupyter` for notebook usage).
  To set up the environment:

  ```bash
  python -m venv .venv
  source .venv/bin/activate
  pip install -r requirements.txt
  ```

  After activation you can run the helper scripts using `python` from the
  same environment.

## Quick start

### 1. Run the tutorial workflow

```bash
python markeR_to_python.py --tutorial --output tutorial.png
open tutorial.png
```

This executes an example from the markeR vignette; consult the
original R tutorial (link below) for a step‑by‑step description of each
analysis step.

### 2. Call any markeR function

```bash
python run_marker_function.py PlotScores \
  --data counts_example \
  --metadata metadata_example \
  --gene_sets genesets_example \
  --Variable "Condition" \
  --method logmedian \
  --nrow 1 \
  --output my_plot.png
```

All parameters are passed as `--name value` pairs.  Use `--help-function` or
refer to the online reference manual for argument names.  The two scripts
are described briefly below.

---

## Tutorial script (`markeR_to_python.py`)

A minimal Python example illustrating how to call markeR from within a Python 
environment. It loads the example datasets, runs a selected function, and 
generates a corresponding plot. The script is intended as a template for users 
who wish to structure and extend their own analysis workflows in Python.

```bash
python markeR_to_python.py --tutorial [--output file.png]
```

---

## Flexible CLI caller (`run_marker_function.py`)

This wrapper constructs and executes R code on the fly, so you can run any
markeR function without _writing_ R.  It handles type conversion, data
loading and optional PNG output.

Basic syntax:

```bash
python run_marker_function.py FUNCTION_NAME [OPTIONS]
```

Useful options:

* `--param value` – named argument for the R function
* `--output filename.png` – capture plot output
* `--width` / `--height` – PNG dimensions (pixels; default 800×600)
* `--help-function` – show documentation link/usage hint
* `--verbose` – display generated R code prior to execution

Example (scores only, no plot):

```bash
python run_marker_function.py CalculateScores \
  --data counts_example \
  --metadata metadata_example \
  --gene_sets genesets_example \
  --method logmedian
```

Example saving a plot to PNG (with custom dimensions):

```bash
python run_marker_function.py PlotScores \
  --data counts_example \
  --metadata metadata_example \
  --gene_sets genesets_example \
  --Variable "Condition" \
  --method logmedian \
  --nrow 1 \
  --width 800 \
  --height 400 \
  --output my_plot.png
```
---

## Built‑in example data

`counts_example`, `metadata_example` and `genesets_example` are loaded
automatically and mirror the objects used in the tutorial.  They let you
try commands without supplying your own datasets.

---

## Tips & troubleshooting

* Add `--verbose` to see the exact R code being run – handy when a
  parameter doesn’t behave as expected.
* If `--output` fails, check write permissions and ensure the directory
  exists.
* Installation problems usually indicate R is missing; make sure `R` is
  on your path before installing Python dependencies.

---

## References

* [markeR on Bioconductor](https://bioconductor.org/packages/markeR)
* [Official R tutorial](https://diseasetranscriptomicslab.github.io/markeR/articles/Article_BenchmarkingMode.html)
* [Reference manual](https://diseasetranscriptomicslab.github.io/markeR/reference/)
* Paper: https://www.biorxiv.org/content/10.64898/2025.12.05.692517

---

*Notes:* scripts install markeR via BiocManager if it’s not already present.
Warnings from `ggplot2` (e.g. about `aes_string()`) are harmless and come
from the package itself.
