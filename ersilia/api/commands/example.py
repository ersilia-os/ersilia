import csv
import json

from ... import ModelBase
from ...core.session import Session
from ...io.input import ExampleGenerator
from ..echo import echo


def example(
    model, simple=True, random=True, n_samples=5, deterministic=False
):
    """
    This command can sample inputs for a given model and save them as a CSV file.

    Args
    -------
    model: The model ID to be served. Can either be the eos identifier or slug.
    file_name: Path where the CSV examples should be saved. Must end with .csv.
    simple: If True, only SMILES strings are returned. If False, outputs include InChiKey and name.
    random: If the model source contains an example input file, when the predefined flag is set, then inputs are sampled from that file. Only the number of samples present in the file are returned, especially if --n_samples is greater than that number. By default, Ersilia samples inputs randomly.
    n_samples: Specify the number of example inputs to generate for the given model.
    deterministic: Used to generate examples data deterministically instead of random sampling. This allows when every time you run with example command with this flag you get the same types of examples.

    Returns
    -------
    Function: The exmaple command function to be used by the API.
    Str: Error message if no model was served in the current session.


    """
    if model:
        model_id = ModelBase(model).model_id
    else:
        session = Session(config_json=None)
        model_id = session.current_model_id()

    if not model_id:
        echo(
            "No model found. Please specify a model or serve a model in the current shell.",
            fg="red",
        )
        return
    
    if not file_name or not file_name.endswith('.csv'):
        echo(
            "Please provide a valid CSV filename ending with .csv",
            fg="red",
            bold=True
        )
        return
    
    eg = ExampleGenerator(model_id=model_id)
    example = eg.example(
        n_samples,
        file_name=None,
        simple=simple,
        try_predefined=not random,
        deterministic=deterministic,
    )
    # normalized = []
    # for item in examples:
    #     if isinstance(item, str):
    #         normalized.append(item)
    #     elif isinstance(item, dict):
    #         smi = item.get('input') or item.get('smiles')
    #         if smi:
    #             normalized.append(smi)
    #     else:
    #         normalized.append(str(item))
    # header = ['input']
    # rows = [[sm] for sm in normalized]

    # try:
    #     with open(file_name, mode='w', newline='', encoding='utf-8') as csvfile:
    #         writer = csv.writer(csvfile)
    #         writer.writerow(header)
    #         writer.writerows(rows)
    #         echo(f":check_mark_button: Examples successfully saved to {file_name}", fg="green", bold=True)
    # except Exception as e:
    #     echo(f"Failed to write examples to CSV: {str(e)}", fg="red", bold=True)
    
    return example
