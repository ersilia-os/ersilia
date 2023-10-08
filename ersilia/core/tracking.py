from datetime import datetime
import pandas as pd

class RunTracker:
    """
    This class will be responsible for tracking model runs. It calculates the desired metadata based on a model's
    inputs, outputs, and other run-specific features, before uploading them to Ersilia's Splunk dashboard.

    NOTE: Currently, the Splunk connection is not set up. For now, we will print tracking results to the console.
    """
    time_start = None
    # def __init__(self):
    #     self.time_start = None

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
        print(self.read_csv(result))

        print("Model metadata:", meta)

        time = datetime.now() - self.time_start
        print("Time taken:", time)
