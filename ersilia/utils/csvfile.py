import os
import csv
import json


class CsvDataLoader(object):
    def __init__(self):
        self.values = None
        self.keys = None
        self.inputs = None
        self.features = None

    def load(self, csv_file):
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
        Reads a file and returns the data as a list of dictionaries.

        :param file_path: Path to the CSV file.
        :return: A list of dictionaries containing the CSV data.
        """

        file_extension = os.path.splitext(file_path)[1].lower()
        if file_extension == ".json":
            return self._read_json(file_path)
        elif file_extension in [".csv", ".tsv"]:
            delimiter = "\t" if file_extension == ".tsv" else ","
            return self._read_csv_tsv(file_path, delimiter)
        else:
            raise ValueError("Unsupported file format")
