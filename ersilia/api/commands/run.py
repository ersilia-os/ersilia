import tempfile

import pandas as pd

from ... import ErsiliaModel
from ...core.session import Session


# input is list
def run(model_id, input, batch_size=100):
    # Runs the current model on a list of SMILES strings and
    # returns the prediction as a pandas data frame.

    #   Args:
    # 		input - a list of SMILES strings
    # 		batch_size - number of SMILES to process per batch

    # 	Returns:
    # 		A pandas df with the predictions
    session = Session(config_json=None)
    model_id = session.current_model_id()
    service_class = session.current_service_class()
    if model_id is None:
        print(
            "No model seems to be served. Please run 'ersilia serve ...' before.",
            fg="red",
        )
        return
    input_df = pd.DataFrame({"input": input})
    with (
        tempfile.NamedTemporaryFile(
            mode="w+t", encoding="utf-8", suffix=".csv", delete=False
        ) as input_file,
        tempfile.NamedTemporaryFile(
            mode="r", encoding="utf-8", suffix=".csv", delete=False
        ) as output_file,
    ):
        mdl = ErsiliaModel(
            model_id,
            output_source=output_file,
            service_class=service_class,
            config_json=None,
        )

        # Question - Formatting of input for run
        input_df.to_csv(input_file.name, index=False)
        input_file.flush()

        result = mdl.run(
            input=input_file.name, output=output_file.name, batch_size=batch_size
        )

        print(
            f"✅ The output successfully generated in {output_file.name} file!",
        )

        output_file.seek(0)
        return pd.read_csv(output_file.name)
