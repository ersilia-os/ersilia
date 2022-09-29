import csv


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
