from datetime import datetime
import json
import pandas as pd


def read_csv(file):
    # reads csv file and returns Pandas dataframe
    return pd.read_csv(file)


class RunTracker:
    """
    This class will be responsible for tracking model runs. It calculates the desired metadata based on a model's
    inputs, outputs, and other run-specific features, before uploading them to Ersilia's Splunk dashboard.

    NOTE: Currently, the Splunk connection is not set up. For now, we will print tracking results to the console.
    """

    def sample_df(self, df, num_rows, num_cols):
        """
        Returns a sample of the dataframe, with the specified number of rows and columns.
        """
        return df.sample(num_rows, axis=0).sample(num_cols, axis=1)

    def __init__(self):
        self.time_start = None

    # function to be called before model is run
    def start_tracking(self):
        self.time_start = datetime.now()

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
        print(read_csv(input))

        print("Run output file:", result)
        resultDf = read_csv(result)
        print(resultDf)

        print("Model metadata:", meta)

        model_id = meta["metadata"].get("Identifier", "Unknown")
        print("Model ID:", model_id)

        time = datetime.now() - self.time_start
        print("Time taken:", time)

        # checking for mismatched types
        nan_count = resultDf.isna().sum()
        print("\nNAN Count:\n", nan_count)

        self.check_types(resultDf, meta["metadata"])

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

    def check_types(self, resultDf, metadata):
        typeDict = {"float64": "Float", "int64": "Int"}
        count = 0

        # ignore key and input columns
        dtypesLst = resultDf.loc[:, ~resultDf.columns.isin(["key", "input"])].dtypes

        for i in dtypesLst:
            if typeDict[str(i)] != metadata["Output Type"][0]:
                count += 1

        if len(dtypesLst) > 1 and metadata["Output Shape"] != "List":
            print("Not right shape. Expected List but got Single")
        elif len(dtypesLst) == 1 and metadata["Output Shape"] != "Single":
            print("Not right shape. Expected Single but got List")
        else:
            print("Output is correct shape.")

        print("Output has", count, "mismatched types.\n")
