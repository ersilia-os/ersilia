"""
Checks of user input shared by the CLI and the Python API.

Each problem is reported with what is wrong and how to fix it. In the CLI the
command then exits with code 1; in the Python API an ``InvalidOptionError``
is raised instead.
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
    Report an error (and a hint): exit with code 1, or raise in the Python API.

    Parameters
    ----------
    message : str
        What went wrong.
    hint : str, optional
        What to do about it.

    Raises
    ------
    InvalidOptionError
        In library mode (the Python API), instead of printing and exiting.
    """
    from .exceptions_utils.throw_ersilia_exception import is_library_mode

    if is_library_mode():
        from .exceptions_utils.api_exceptions import InvalidOptionError

        raise InvalidOptionError(message, hint or "")
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
    _check_input_file(input)
    if output is None or not output.lower().endswith((".csv", ".h5")):
        fail("The output file must end in .csv or .h5.")
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


def check_fetch_folder(model, from_dir):
    """
    Check a model folder given to fetch (``--from_dir``).

    Parameters
    ----------
    model : str
        The model the user asked for (identifier or slug).
    from_dir : str
        The folder with the model repository.

    Returns
    -------
    ModelBase
        The model found in the folder.
    """
    from ..core.modelbase import ModelBase
    from .paths import get_metadata_from_base_dir

    if not os.path.isdir(os.path.expanduser(from_dir)):
        fail(f"The folder {from_dir} does not exist.")
    mdl = ModelBase(repo_path=from_dir)
    try:
        folder_id = get_metadata_from_base_dir(from_dir).get("Identifier")
    except Exception:
        folder_id = None
    folder_id = folder_id or mdl.model_id
    if folder_id and model.strip().lower() not in (folder_id, mdl.slug):
        fail(
            f"The folder {from_dir} contains model {folder_id}, not {model}.",
            f"Run 'ersilia fetch {folder_id} --from_dir {from_dir}'.",
        )
    return mdl


def check_curated_examples(model_id):
    """
    Check that a fetched model has its own (curated) example inputs.

    Parameters
    ----------
    model_id : str
        The model identifier.
    """
    from ..core.modelbase import ModelBase
    from ..default import PREDEFINED_EXAMPLE_FILES

    model_dir = ModelBase(model_id)._model_path(model_id)
    if not any(
        os.path.exists(os.path.join(model_dir, f)) for f in PREDEFINED_EXAMPLE_FILES
    ):
        fail(
            f"Model {model_id} has no curated examples here.",
            "Fetch the model first, or use the random mode instead.",
        )
