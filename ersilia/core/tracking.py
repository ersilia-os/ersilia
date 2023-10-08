class RunTracker:
    """
    This class will be responsible for tracking model runs. It calculates the desired metadata based on a model's
    inputs, outputs, and other run-specific features, before uploading them to Ersilia's Splunk dashboard.

    NOTE: Currently, the Splunk connection is not set up. For now, we will print tracking results to the console.
    """

    def track(self, input, result, meta):
        """
        Tracks the results after a model run.
        """

        print("Run input file:", input)
        print("Run output file:", result)

        print("Model metadata:", meta)
