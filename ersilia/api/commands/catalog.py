import io
import json

import pandas as pd

from ... import logger
from ...hub.content.card import ModelCard
from ...hub.content.catalog import ModelCatalog
from ..echo import echo


def catalog(
    hub=False,
    file_name=None,
    more=False,
    card=False,
    model=None,
    as_json=False,
    verbose=False,
):
    """
    API-compatible version of the catalog command with echo-based output.

    This command allows users to list a catalog of models available either locally or in the model hub.
    It provides options to display the catalog in various formats(such as tables by default or json), show more detailed information,
    and view model cards for specific models.

    Parameters
    ----------
    hub : bool
        If True, fetch the catalog from the hub.
    file_name : str or None
        If specified, write the catalog to this file.
    more : bool
        If True, show more detail in catalog.
    card : bool
        If True, display the model card for a given model.
    model : str or None
        The model ID for which to display metadata.
    as_json : bool
        If True, return JSON output instead of a formatted table.
    verbose : bool
        If True, enable verbose logging.

    Returns
    -------
    pandas.DataFrame or dict or None
        A DataFrame of the last two catalog columns (for non-card mode),
        the model card metadata as dict (for card mode),
        or None on failure.
    """
    logger.set_verbosity(1 if verbose else 0)

    if card:
        if not model:
            echo("❌ Error: --card option requires a model ID", fg="red", bold=True)
            # return None
        try:
            mc = ModelCard()
            metadata = mc.get(model, as_json=True)

            if not metadata:
                echo(
                    f"❌ Error: No metadata found for model ID '{model}'",
                    fg="red",
                    bold=True,
                )
                # return None

            if as_json:
                echo(json.dumps(metadata, indent=2))
            else:
                echo(str(metadata))
            return metadata

        except Exception as e:
            echo(f"❌ Error fetching model metadata: {e}", fg="red", bold=True)
            # return None

    # Not in card mode → fetch catalog
    catalog_obj = ModelCatalog(less=not more)
    try:
        catalog_table = catalog_obj.hub() if hub else catalog_obj.local()
    except Exception as e:
        echo(f"❌ Error loading catalog: {e}", fg="red", bold=True)

    if not catalog_table.data:
        echo(
            "⚠️ No local models found. Try `ersilia fetch` to download a model.",
            fg="yellow",
        )

    # Save to file if requested (but do NOT return early)
    if file_name:
        catalog_table.write(file_name)
        echo(f"📁 Catalog written to {file_name}", fg="green")

    df = pd.read_json(io.StringIO(catalog_table.as_json()))
    if df.shape[0] > 0:
        df = df.drop(columns=["Index"])
        return df
    else:
        return None
