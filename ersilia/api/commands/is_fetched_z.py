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
        return False

    if not catalog_table.data:
        echo("❌No models are fetched locally", fg="red", bold=True)
        return False

    # Check if "Identifier" column exists in the catalog table
    if "Identifier" not in catalog_table.columns:
        echo(
            "❌Unexpected catalog structure: 'Identifier' column not found",
            fg="red",
            bold=True,
        )
        return False

    # Get the index of the Identifier column
    identifier_idx = catalog_table.columns.index("Identifier")

    # Check if the model_id exists in the Identifier column
    for row in catalog_table.data:
        if row[identifier_idx] == model_id:
            echo(f"✅ Model {model_id} is fetched", fg="green", bold=True)
            return True

    echo(f"❌Model {model_id} is not fetched", fg="red", bold=True)
    return False
