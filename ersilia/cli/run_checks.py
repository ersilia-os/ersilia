"""
Checks of `ersilia run` inputs and outputs, done before the model is loaded.

Each problem is reported with what is wrong and how to fix it, and the command
exits with code 1.
"""

import csv
import os
import re
import sys
import tempfile

from .echo import echo

INPUT_HEADERS = {"smiles", "input"}
# A first line made only of SMILES atoms and symbols (e.g. "CCO", "c1ccccc1")
# is an input, not a header. Ordinary column names ("name", "smiles", "ID")
# do not match.
_LOOKS_LIKE_SMILES = re.compile(r"[CNOSPFIBrclnosp0-9()=#\[\]@+\-\\/.%]{2,}")


def fail(message, hint=None):
    """
    Print an error (and a hint) and exit with code 1.

    Parameters
    ----------
    message : str
        What went wrong.
    hint : str, optional
        What to do about it.
    """
    echo(message, fg="red")
    if hint:
        echo(hint)
    sys.exit(1)


def _check_input_file(path):
    if path.lower().endswith((".xlsx", ".xls")):
        fail(
            f"{path} is an Excel file; the input must be a CSV file.",
            "Export it as CSV (one column of inputs) first.",
        )
    if not path.lower().endswith(".csv"):
        if not os.path.exists(path) and os.sep not in path and "." not in path[-5:]:
            fail(
                "The input must be a CSV file, not a single input.",
                "To run one molecule, put it in a CSV file, e.g. printf 'smiles\\nCCO\\n' > input.csv",
            )
        fail(
            "The input must be a CSV file with one column of inputs.",
            "Save it as a comma-separated .csv file.",
        )
    if not os.path.exists(path):
        fail(f"Input file {path} does not exist.")
    if os.path.isdir(path):
        fail(f"{path} is a folder, not a file.")
    if os.path.getsize(path) == 0:
        fail(
            f"The input file {path} is empty.",
            "It needs a header line (e.g. 'smiles') and one input per line.",
        )


def _read_rows(path):
    try:
        with open(path, "r", encoding="utf-8-sig", newline="") as f:
            return list(csv.reader(f))
    except UnicodeDecodeError:
        fail(
            f"The input file {path} is not UTF-8 encoded.",
            'Save it as "CSV UTF-8" and try again.',
        )
    except csv.Error as e:
        fail(f"The input file {path} could not be read as CSV: {e}.")


def check_run_arguments(input, output):
    """
    Check the input and output of `ersilia run`, and prepare the input.

    Parameters
    ----------
    input : str
        Path of the input CSV file.
    output : str
        Path of the output file (.csv or .h5).

    Returns
    -------
    tuple of (str, str or None)
        The input path to run, and a temporary folder to remove afterwards
        (when a single column was extracted from a wider file), or None.
    """
    from .messages import wrong_extension

    _check_input_file(input)
    if output is None or not output.lower().endswith((".csv", ".h5")):
        wrong_extension([".csv", ".h5"])
        sys.exit(1)
    if os.path.realpath(input) == os.path.realpath(output):
        fail(
            "The output file is the same as the input file.",
            "Choose a different output file, so the input is not overwritten.",
        )
    folder = os.path.dirname(output) or "."
    if not os.path.isdir(folder):
        fail(f"The output folder {folder} does not exist.")
    if not os.access(folder, os.W_OK):
        fail(f"Cannot write to the output folder {folder} (permission denied).")

    rows = [r for r in _read_rows(input) if any(c.strip() for c in r)]
    if not rows:
        fail(
            f"The input file {input} has no content.",
            "It needs a header line (e.g. 'smiles') and one input per line.",
        )
    header, data = rows[0], rows[1:]
    if len(header) == 1 and ("\t" in header[0] or ";" in header[0]):
        fail(
            "The input file looks tab or semicolon separated.",
            "Save it as a comma-separated CSV file with one column of inputs.",
        )
    if not data:
        fail(
            f"The input file {input} has a header but no inputs.",
            "Add one input per line below the header.",
        )

    tmp_dir = None
    if len(header) > 1:
        named = [i for i, h in enumerate(header) if h.strip().lower() in INPUT_HEADERS]
        columns = ", ".join(h.strip() for h in header)
        if len(named) != 1:
            fail(
                f"The input file has {len(header)} columns ({columns}).",
                "Keep only the column with the inputs, or name it 'smiles'.",
            )
        col = named[0]
        echo(
            f"The input file has {len(header)} columns; using '{header[col].strip()}'.",
            fg="yellow",
        )
        tmp_dir = tempfile.mkdtemp(prefix="ersilia-run-")
        single = os.path.join(tmp_dir, "input.csv")
        with open(single, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([header[col]])
            for row in data:
                writer.writerow([row[col] if col < len(row) else ""])
        input = single
    else:
        first = header[0].strip()
        if first.lower() not in INPUT_HEADERS and _LOOKS_LIKE_SMILES.fullmatch(first):
            echo(
                f"The first line '{first}' looks like an input, not a header. It was treated as a header and skipped.",
                fg="yellow",
            )
            echo("Add a header line (e.g. 'smiles') to include it.")
    return input, tmp_dir
