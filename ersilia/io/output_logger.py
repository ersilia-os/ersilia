import csv
import os


class TabularResultLogger(object):
    """
    Class to log tabular results.

    This class handles logging of results from tabular files such as CSV and TSV.

    Parameters
    ----------
    None
    """

    def __init__(self):
        pass

    def _is_tabular_file(self, s):
        if type(s) is str:
            if not os.path.exists(s):
                return False
            if s.endswith(".csv"):
                return True
            if s.endswith(".tsv"):
                return True
            if s.endswith(".h5"):
                return True
        else:
            return False

    def _get_delimiter(self, s):
        if s.endswith(".csv"):
            return ","
        if s.endswith(".tsv"):
            return "\t"
        return None

    def summary(self, output):
        """
        Summarize results statistics
        This class produces a summary of the results, counting inputs, output dimensions, etc.
        Parameters
        ----------
        output : str
            Path to the output file.
        Returns
        -------
        data : dict
            Summary of the results.
        """
        if not self._is_tabular_file(output):
            raise ValueError("Output file is not a valid tabular file.")

        if output.endswith(".h5"):
            raise ValueError("Summarization is not yet available for HDF5 files.")

        delimiter = self._get_delimiter(output)

        with open(output, "r") as f:
            reader = csv.reader(f, delimiter=delimiter)
            output_dim = len(next(reader)) - 2
            num_inputs = 0
            full_empty_rows = 0
            num_empty_cells = 0
            for row in reader:
                row = row[2:]
                num_empty_cells_ = sum(
                    1
                    for cell in row
                    if cell == "" or cell == "nan" or cell == "NaN" or cell == "None"
                )
                if num_empty_cells_ == output_dim:
                    full_empty_rows += 1
                num_empty_cells += num_empty_cells_
                num_inputs += 1
            proportion_empty_cells = num_empty_cells / (num_inputs * output_dim)

        data = {
            "num_inputs": num_inputs,
            "output_dim": output_dim,
            "prop_empty_cells": proportion_empty_cells,
            "full_empty_rows": full_empty_rows,
        }

        return data
