import json

from ...hub.content.card import ModelCard
from ...hub.content.catalog import ModelCatalog
from ... import ErsiliaModel, logger
from ..echo import echo
import pandas as pd



# def catalog_cmd():
#     """
#     Creates the catalog command for the API.

#     This command allows users to list a catalog of models available either locally or in the model hub.
#     It provides options to display the catalog in various formats(such as tables by default or json), show more detailed information,
#     and view model cards for specific models.

#     Returns
#     -------
#     function
#         The catalog command function to be used by the API and for testing in the pytest.

#     """

def catalog(
        hub=False,
        file_name=None,
        browser=False,
        more=False,
        card=False,
        model=None,
        as_json=False,
        verbose=False,
):
        """
        API-compatible version of the catalog command with echo-based output.

        Parameters
        ----------
        hub : bool
            If True, fetch the catalog from the hub.
        file_name : str or None
            If specified, write the catalog to this file.
        browser : bool
            Unused in current version; reserved for future.
        more : bool
            If True, show more detail in catalog.
        card : bool
            If True, display the model card for a given model.
        model : str or None
            The model ID for which to display metadata.
        as_json : bool
            If True, return JSON output instead of a formatted table.

        Returns
        -------
        pandas.DataFrame
        A DataFrame containing the last two columns of the model catalog.
        Also prints the full catalog as a table or JSON to the terminal, depending on `as_json`.
        """
        if verbose:
            logger.set_verbosity(1)
        else:
            logger.set_verbosity(0)

        if card and not model:
            echo(" Error: --card option requires a model ID", err=True)
            return None

        if card and model:
            try:
                mc = ModelCard()
                model_metadata = mc.get(model, as_json=True)

                if not model_metadata:
                    echo(f" Error: No metadata found for model ID '{model}'", err=True)
                    return None

                if as_json:
                    echo(json.dumps(model_metadata, indent=2))
                    return model_metadata
                else:
                    echo(str(model_metadata))
                    return model_metadata

            except Exception as e:
                echo(f"Error fetching model metadata: {e}", err=True)
                return None

        else:
            mc = ModelCatalog(less=not more)

            if hub:
                catalog_table = mc.hub()
            else:
                catalog_table = mc.local()
                if not catalog_table.data:
                    echo(
                        " No local models available. Please fetch a model by running 'ersilia fetch'",
                        err=True,
                    )
                    return None
            
            if file_name is None:
                catalog = (
                    catalog_table.as_json() if as_json else catalog_table.as_table()
                )

            else:
                catalog_table.write(file_name)
                echo(f"📁 Catalog written to {file_name}")
                return None
            
            try:
                df = pd.DataFrame(catalog_table.data)
                if df.shape[1] >= 2:
                    df = df.iloc[:, -2:] 
            except Exception as e:
                echo(f"❌ Could not convert catalog to DataFrame: {e}", err=True)

        echo(catalog)
        return df
