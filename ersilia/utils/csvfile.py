import csv
import json
import os


class CsvDataLoader(object):
    """
    A class to load and process data from CSV and JSON files.

    Methods
    -------
    load(csv_file)
        Load data from a CSV file.
    read(file_path)
        Read data from a CSV, TSV, or JSON file.
    """

    def __init__(self):
        self.values = None
        self.keys = None
        self.inputs = None
        self.features = None

    def load(self, csv_file):
        """
        Load data from a CSV file.

        Parameters
        ----------
        csv_file : str
            The path to the CSV file.
        """
        with open(csv_file, "r") as f:
            reader = csv.reader(f)
            self.features = [
                x for x in next(reader) if x not in ["key", "input", "text"]
            ]
            self.keys = []
            self.inputs = []
            self.values = []
            for r in reader:
                self.keys += [r[0]]
                self.inputs += [r[1]]
                self.values += [r[-len(self.features) :]]

    def _read_csv_tsv(self, file_path, delimiter):
        with open(file_path, mode="r") as file:
            reader = csv.DictReader(file, delimiter=delimiter)
            data = [row for row in reader]
            return data

    def _read_json(self, file_path):
        with open(file_path, mode="r") as file:
            return json.load(file)

    def read(self, file_path):
        """
        Read data from a CSV, TSV, or JSON file.

        Parameters
        ----------
        file_path : str
            The path to the file.

        Returns
        -------
        list or dict
            The data read from the file.

        Raises
        ------
        ValueError
            If the file format is unsupported.
        """
        file_extension = os.path.splitext(file_path)[1].lower()
        if file_extension == ".json":
            return self._read_json(file_path)
        elif file_extension in [".csv", ".tsv"]:
            delimiter = "\t" if file_extension == ".tsv" else ","
            return self._read_csv_tsv(file_path, delimiter)
        else:
            raise ValueError("Unsupported file format")
