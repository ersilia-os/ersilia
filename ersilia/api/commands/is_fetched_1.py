import io

import pandas as pd

from ...hub.content.catalog import ModelCatalog
from ..echo import echo


def is_fetched(model_id: str) -> bool:
    """
    Check whether a given Ersilia model ID has been fetched locally.

    Parameters
    ----------
    model_id : str
        The identifier of the model to check
    Returns
    ------
    bool
        True if 'model_id' exists in the local catalog (i.e., has been fetched), False otherwise.
    """
    catalog_obj = ModelCatalog(less=True)
    try:
        catalog_table = catalog_obj.local()
    except Exception as e:
        echo(f"❌Error loading local catalog: {e}", fg="red", bold=True)

    if not catalog_table.data:
        return False
    df = pd.read_json(io.StringIO(catalog_table.as_json()))
    if "Identifier" in df.columns:
        return model_id in df["Identifier"].values
    return False
