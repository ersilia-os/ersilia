from datetime import datetime
import json
import pandas as pd


class RunTracker:
    """
    This class will be responsible for tracking model runs. It calculates the desired metadata based on a model's
    inputs, outputs, and other run-specific features, before uploading them to Ersilia's Splunk dashboard.

    NOTE: Currently, the Splunk connection is not set up. For now, we will print tracking results to the console.
    """

    def __init__(self):
        self.time_start = None

    # function to be called before model is run
    def start_tracking(self):
        self.time_start = datetime.now()

    def read_csv(self, file):
        # reads csv file and returns Pandas dataframe
        return pd.read_csv(file)

    def stats(self, result):
        dat = self.read_csv(result)

        # drop first two columns (key, input)
        dat = dat.drop(["key", "input"], axis=1)

        # calculate and print statistics
        for column in dat:
            print("Mean %s: %s" % (column, dat[column].mean()))
            if len(dat[column].mode()) == 1:
                print("Mode %s: %s" % (column, dat[column].mode()))
            else:
                print("No mode")
            print("Min %s: %s" % (column, dat[column].min()))
            print("Max %s: %s" % (column, dat[column].max()))
            print("Standard deviation %s: %s" % (column, dat[column].std()))

    def track(self, input, result, meta):
        """
        Tracks the results after a model run.
        """
        print("Run input file:", input)
        print(self.read_csv(input))

        print("Run output file:", result)
        print(self.read_csv(result))

        print("Model metadata:", meta)

        model_id = meta["metadata"].get("Identifier", "Unknown")
        print("Model ID:", model_id)

        time = datetime.now() - self.time_start
        print("Time taken:", time)

        self.stats(result)

        input_dataframe = self.read_csv(input)
        result_dataframe = self.read_csv(result)

        input_size = input_dataframe.memory_usage(deep=True).sum() / 1024
        output_size = result_dataframe.memory_usage(deep=True).sum() / 1024

        input_avg_row_size = input_size / len(input_dataframe)
        output_avg_row_size = output_size / len(result_dataframe)

        print("Average Input Row Size (KB):", input_avg_row_size)
        print("Average Output Row Size (KB):", output_avg_row_size)

    def log_to_console(self, data):
        print(f"\n{json.dumps(data)}\n")

    def read_json(self, result):
        data = json.load(result)
        self.log_to_console(result)
        return data
