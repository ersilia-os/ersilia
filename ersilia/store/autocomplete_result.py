import csv
import json
import os
import tempfile

from ..utils.terminal import run_command


class CalculateMissingResult:
    """
    A class to identify missing values in a CSV, compute replacements via an external service,
    and update the CSV efficiently.

    Parameters
    ----------
    path : str
        Path to the input CSV file.
    temp_output_path : str
        Path to temporary output file for computation results.
    model_id : str
        Identifier of the external computation model.
    """

    def __init__(self, path, temp_output_path, model_id, cache_saving):
        self.cache_saving = cache_saving
        self.model_id = model_id
        self.path = path
        self.temp_output_path = temp_output_path
        self.missing_rows = []
        self.results = []
        self.header = []

    def find_missing(self):
        """
        Identify rows where the third column is empty.

        Returns
        -------
        list of dict
            Each dict contains 'key', 'input', and 'index' of a missing row.
        """
        with open(self.path, newline="", encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile)
            self.header = next(reader)
            for idx, row in enumerate(reader, start=1):
                if len(row) < 3 or not row[2].strip():
                    self.missing_rows.append(
                        {"key": row[0], "input": row[1], "index": idx}
                    )
        return self.missing_rows

    def compute(self):
        """
        Compute missing values by batching inputs and invoking the service.

        Returns
        -------
        list of list
            Lists of computed values corresponding to missing_rows order.
        """
        inputs = [item["input"] for item in self.missing_rows]
        inputs_json = json.dumps(inputs)
        self._run(inputs_json)
        self.results = self._read_temp_result()
        return self.results

    def update_csv(self, output_csv_path=None):
        """
        Write a new CSV with computed values inserted at the third column onward.

        Parameters
        ----------
        output_csv_path : str, optional
            Path for the updated CSV file. Overwrites original if None or same as source.

        Returns
        -------
        str
            Path to the updated CSV file.
        """
        src = self.path
        dest_given = (
            output_csv_path if output_csv_path and output_csv_path != src else None
        )
        dst = dest_given or f"{src}.tmp"

        with (
            open(src, newline="", encoding="utf-8") as infile,
            open(dst, "w", newline="", encoding="utf-8") as outfile,
        ):
            reader = csv.reader(infile)
            writer = csv.writer(outfile)
            writer.writerow(self.header)

            next(reader, None)
            for idx, row in enumerate(reader, start=1):
                if idx in {info["index"] for info in self.missing_rows}:
                    outputs = {
                        info["index"]: vals
                        for info, vals in zip(self.missing_rows, self.results)
                    }[idx]
                    needed = 2 + len(outputs) - len(row)
                    if needed > 0:
                        row.extend([""] * needed)
                    for i, val in enumerate(outputs):
                        row[2 + i] = val
                writer.writerow(row)

        if dest_given is None:
            os.replace(dst, src)
        os.remove(self.temp_output_path)
        return output_csv_path or src

    def _run(self, inputs_json):
        inputs = self._make_smiles_csv(inputs_json)
        cmd = (
            f"ersilia -v serve {self.model_id} --cache-saving {self.cache_saving} && "
            f"ersilia -v run -i '{inputs}' -o '{self.temp_output_path}'"
        )
        return run_command(cmd, quiet=True)

    def _read_temp_result(self):
        try:
            with open(self.temp_output_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            return [list(item) for item in data]
        except (json.JSONDecodeError, ValueError):
            with open(self.temp_output_path, "r", encoding="utf-8") as f:
                reader = csv.reader(f)
                next(reader, None)
                return [row[2:] for row in reader]

    def _make_smiles_csv(self, smiles_list):
        smiles_list = json.loads(smiles_list)
        tmp = tempfile.NamedTemporaryFile(
            mode="w", newline="", suffix=".csv", delete=False
        )
        writer = csv.writer(tmp)
        writer.writerow(["input"])
        for smi in smiles_list:
            writer.writerow([smi])

        tmp.flush()
        tmp.close()
        return tmp.name

    def process(self, output_csv_path=None):
        """
        Full pipeline: find missing, compute, and update CSV.

        Parameters
        ----------
        output_csv_path : str, optional
            Path for the updated CSV file.

        Returns
        -------
        str
            Path to the updated CSV file.
        """
        self.find_missing()
        self.compute()
        return self.update_csv(output_csv_path)
