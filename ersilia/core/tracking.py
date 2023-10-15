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

    def track(self, input, result, meta):
        """
        Tracks the results after a model run.
        """

        print("Run input file:", input)
        print(self.read_csv(input))

        print("Run output file:", result)
        resultDf = self.read_csv(result)
        print(resultDf)

        print("Model metadata:", meta)

        time = datetime.now() - self.time_start
        print("Time taken:", time)

        # checking for mismatched types
        nan_count = resultDf.isna().sum()
        print("\nNAN Count:\n", nan_count)

        self.check_types(resultDf, meta['metadata'])

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
        dtypesLst = resultDf.loc[:, ~resultDf.columns.isin(['key', 'input'])].dtypes

        for i in dtypesLst:
            if typeDict[str(i)] != metadata['Output Type'][0]:
                count += 1

        if (len(dtypesLst) > 1 and metadata['Output Shape'] != 'List'):
            print("Not right shape. Expected List but got Single")
        elif (len(dtypesLst) == 1 and metadata['Output Shape'] != 'Single'):
            print("Not right shape. Expected Single but got List")
        else:
            print("Output is correct shape.")

        print("Output has", count, "mismatched types.\n")

        

