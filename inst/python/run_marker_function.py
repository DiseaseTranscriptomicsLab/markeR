"""Flexible wrapper to call any markeR R function from Python.

This script allows you to call markeR functions directly without writing Python code.

Usage examples:
    python run_marker_function.py PlotScores --help-function
    python run_marker_function.py CalculateScores \\
        --data counts_example --metadata metadata_example \\
        --gene_sets genesets_example --method logmedian --verbose
    
    python run_marker_function.py PlotScores \\
        --data counts_example --metadata metadata_example \\
        --gene_sets genesets_example --Variable "Condition" \\
        --method logmedian --nrow 1 --output plot.png

Options:
    --help-function      Show R documentation for the function
    --verbose            Print the generated R code before executing
    --output FILE        Save plot to PNG file

For built-in example data, use the names: counts_example, metadata_example, genesets_example
"""

import sys
import argparse
import json
import os

# Check dependencies
_missing = []
try:
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr, isinstalled
except ImportError:
    _missing.append("rpy2")

if _missing:
    sys.exit(
        "The following Python packages are required but not installed: %s.\n"
        "Please install them (e.g. `pip install rpy2`)." % ", ".join(_missing)
    )

def ensure_bioc_installed() -> None:
    """Install Bioconductor's package manager if it is not already present."""
    utils = importr("utils")
    biocinstaller = "BiocManager"
    if not isinstalled(biocinstaller):
        ro.r('install.packages("{0}")'.format(biocinstaller))
    ro.r('suppressMessages(require({0}))'.format(biocinstaller))


def install_markeR() -> None:
    """Install the markeR package from Bioconductor if not already installed."""
    ensure_bioc_installed()
    try:
        importr("markeR")
    except Exception:
        ro.r('BiocManager::install("markeR", ask=FALSE, update=FALSE)')
        ro.r('library(markeR)')


def load_example_data():
    """Load built-in markeR example datasets into R namespace."""
    install_markeR()
    # Load the example datasets
    ro.r('data("genesets_example", package="markeR")')
    ro.r('data("counts_example", package="markeR")')
    ro.r('data("metadata_example", package="markeR")')
    print("Loaded markeR example datasets: counts_example, metadata_example, genesets_example")


def parse_parameter(value: str):
    """
    Parse a parameter value intelligently.
    - Numbers become numeric
    - "true"/"false" become logical
    - "null" becomes NULL
    - R object names (e.g., counts_example) are kept as-is
    - JSON objects/arrays become R equivalents
    - Strings are kept as strings
    """
    value_lower = value.lower()
    
    # Handle boolean
    if value_lower == "true":
        return "TRUE"
    if value_lower == "false":
        return "FALSE"
    if value_lower == "null":
        return "NULL"
    
    # Handle numbers
    try:
        if "." in value:
            float(value)
            return value
        else:
            int(value)
            return value
    except ValueError:
        pass
    
    # Check if it's a known R object name (example data)
    known_objects = ["counts_example", "metadata_example", "genesets_example"]
    if value in known_objects:
        return value
    
    # Try JSON parsing for objects/arrays
    try:
        json.loads(value)
        # If it parses as JSON, return as-is (user can provide lists as JSON)
        return value
    except (json.JSONDecodeError, ValueError):
        pass
    
    # Default: treat as string
    return f'"{value}"'


def build_r_call(function_name: str, params: dict, output_file: str = None, width: int = 800, height: int = 600) -> str:
    """
    Build an R function call string from parameters.
    
    Parameters
    ----------
    function_name : str
        Name of the R function to call
    params : dict
        Dictionary of parameter names and values
    output_file : str
        If provided, set up PNG device before the call and close after
    width : int
        PNG width in pixels (default: 800)
    height : int
        PNG height in pixels (default: 600)
        
    Returns
    -------
    str
        Complete R code to execute
    """
    # Remove output_file from params if present
    params = {k: v for k, v in params.items() if k != "output_file"}
    
    # Build parameter list
    param_strings = []
    for key, value in params.items():
        parsed_value = parse_parameter(value)
        param_strings.append(f"{key} = {parsed_value}")
    
    param_str = ", ".join(param_strings)
    
    # Build R code
    r_code = ""
    
    if output_file:
        output_file = os.path.abspath(output_file)
        dirname = os.path.dirname(output_file)
        if dirname and not os.path.isdir(dirname):
            os.makedirs(dirname, exist_ok=True)
        r_code += f'png("{output_file}", width={width}, height={height})\n'
        r_code += f"result <- {function_name}({param_str})\n"
        r_code += "tryCatch(print(result), error=function(e) { invisible(NULL) })\n"
        r_code += "dev.off()\n"
    else:
        # For screen output, try to print the result
        r_code += f"result <- {function_name}({param_str})\n"
        r_code += "tryCatch(print(result), error=function(e) { cat('Function executed.\\n') })\n"
    
    return r_code


def main():
    parser = argparse.ArgumentParser(
        description="Call any markeR R function from Python",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_marker_function.py CalculateScores \\
    --data counts_example --metadata metadata_example \\
    --gene_sets genesets_example --method logmedian
    
  python run_marker_function.py PlotScores \\
    --data counts_example --metadata metadata_example \\
    --gene_sets genesets_example --Variable "Condition" \\
    --method logmedian --output my_plot.png

Built-in example data: counts_example, metadata_example, genesets_example
        """
    )
    
    parser.add_argument(
        "function_name",
        help="Name of the markeR function to call (e.g., CalculateScores, PlotScores)"
    )
    
    parser.add_argument(
        "--help-function",
        action="store_true",
        help="Show help for the R function (instead of calling it)"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print the generated R code before executing"
    )
    
    parser.add_argument(
        "--output",
        help="Save plot output to a PNG file"
    )
    
    parser.add_argument(
        "--width",
        type=int,
        default=800,
        help="PNG width in pixels (default: 800)"
    )
    
    parser.add_argument(
        "--height",
        type=int,
        default=600,
        help="PNG height in pixels (default: 600)"
    )
    
    # Allow arbitrary parameters
    parser.add_argument(
        "params",
        nargs="*",
        help="Parameters as --name value pairs (e.g., --data counts_example --method logmedian)"
    )
    
    # Handle --help for specific functions
    if len(sys.argv) > 1 and sys.argv[1] not in ["--help", "-h"]:
        if "--help-function" in sys.argv:
            func_name = sys.argv[1]
            print(f"\n{'='*70}")
            print(f"Help for markeR::{func_name}")
            print(f"{'='*70}\n")
            install_markeR()
            # Try to display help
            try:
                # Get function signature and description
                ro.r(f'''
library(markeR)
cat("Function: {func_name}\\n\\n")
# Try to get help
tryCatch({{
  help_file <- help("{func_name}", package="markeR")
  # Get description from help
}}, error = function(e) {{
  cat("Help available at: https://diseasetranscriptomicslab.github.io/markeR/reference/{func_name}.html\\n")
}})
''')
            except Exception as e:
                pass
            
            print(f"\nDocumentation:")
            print(f"  https://diseasetranscriptomicslab.github.io/markeR/reference/{func_name}.html")
            print(f"\nTo use this function:")
            print(f"  python run_marker_function.py {func_name} --param1 value1 --param2 value2 [--output output.png]")
            print(f"\nTip: Use --verbose flag to see generated R code")
            print(f"  python run_marker_function.py {func_name} --verbose --param1 value1 ...\n")
            print(f"{'='*70}\n")
            return
    
    # Parse args
    if len(sys.argv) < 2:
        parser.print_help()
        return
    
    func_name = sys.argv[1]
    output_file = None
    width = 800
    height = 600
    
    # Parse remaining arguments as key-value pairs
    params = {}
    i = 2
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == "--output" and i + 1 < len(sys.argv):
            output_file = sys.argv[i + 1]
            i += 2
        elif arg == "--width" and i + 1 < len(sys.argv):
            try:
                width = int(sys.argv[i + 1])
            except ValueError:
                print(f"Error: --width must be a number, got '{sys.argv[i + 1]}'")
                sys.exit(1)
            i += 2
        elif arg == "--height" and i + 1 < len(sys.argv):
            try:
                height = int(sys.argv[i + 1])
            except ValueError:
                print(f"Error: --height must be a number, got '{sys.argv[i + 1]}'")
                sys.exit(1)
            i += 2
        elif arg in ["--verbose", "--help-function"]:
            # Skip flags that are not parameters
            i += 1
        elif arg.startswith("--"):
            key = arg[2:]  # Remove --
            if i + 1 < len(sys.argv) and not sys.argv[i + 1].startswith("--"):
                value = sys.argv[i + 1]
                params[key] = value
                i += 2
            else:
                # Boolean flag
                params[key] = "TRUE"
                i += 1
        else:
            i += 1
    
    # Ensure markeR is installed
    print("Installing markeR if needed...")
    install_markeR()
    load_example_data()
    
    # Build and execute the R call
    print(f"\nCalling {func_name} with parameters:")
    for key, value in params.items():
        print(f"  {key} = {value}")
    
    if output_file:
        print(f"  Saving plot to: {output_file}")
        print(f"  PNG dimensions: {width}x{height} pixels")
    
    r_code = build_r_call(func_name, params, output_file, width, height)
    
    # Show R code if verbose mode
    if "--verbose" in sys.argv:
        print(f"\n{'='*70}")
        print("Generated R code:")
        print(f"{'='*70}")
        print(r_code)
        print(f"{'='*70}\n")
    
    print(f"Executing R code...\n")
    print("=" * 60)
    
    try:
        ro.r(r_code)
        print("=" * 60)
        if output_file:
            print(f"\n✓ Plot saved to: {output_file}")
        else:
            print(f"\n✓ Function executed successfully")
    except Exception as e:
        print("=" * 60)
        print(f"\n✗ Error executing function: {e}")
        if "--verbose" not in sys.argv:
            print("\nTip: Use --verbose flag to see the generated R code")
            print(f"  python run_marker_function.py {func_name} --verbose [other options]")
        sys.exit(1)


if __name__ == "__main__":
    main()
