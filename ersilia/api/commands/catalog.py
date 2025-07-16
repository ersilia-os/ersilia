import json

from ...hub.content.card import ModelCard
from ...hub.content.catalog import ModelCatalog
from ... import ErsiliaModel, logger


def catalog_cmd():
    """
    Creates the catalog command for the API.

    This command allows users to list a catalog of models available either locally or in the model hub.
    It provides options to display the catalog in various formats(such as tables by default or json), show more detailed information,
    and view model cards for specific models.

    Returns
    -------
    function
        The catalog command function to be used by the API and for testing in the pytest.

    """

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
        str or dict or None
            Catalog data or model metadata, or None if written to file or silenced.
        """
        if verbose:
            logger.set_verbosity(1)
        else:
            logger.set_verbosity(0)

        if card and not model:
            click.echo("❌ Error: --card option requires a model ID", err=True)
            return None

        if card and model:
            try:
                mc = ModelCard()
                model_metadata = mc.get(model, as_json=True)

                if not model_metadata:
                    click.echo(f"❌ Error: No metadata found for model ID '{model}'", err=True)
                    return None

                if as_json:
                    click.echo(json.dumps(model_metadata, indent=2))
                    return model_metadata
                else:
                    click.echo(str(model_metadata))
                    return model_metadata

            except Exception as e:
                click.echo(f"❌ Error fetching model metadata: {e}", err=True)
                return None

        else:
            mc = ModelCatalog(less=not more)

            if hub:
                catalog_table = mc.hub()
            else:
                catalog_table = mc.local()
                if not catalog_table.data:
                    click.echo(
                        "⚠️ No local models available. Please fetch a model by running 'ersilia fetch'",
                        err=True,
                    )
                    return None

            if file_name is None:
                catalog = (
                    catalog_table.as_json() if as_json else catalog_table.as_table()
                )
                click.echo(catalog)
                return catalog
            else:
                catalog_table.write(file_name)
                click.echo(f"📁 Catalog written to {file_name}")
                return None

    return catalog

    def catalog(
        hub=False,
        file_name=None,
        browser=False,
        more=False,
        card=False,
        model=None,
        as_json=False,
    ):
        if card and not model:
            click.echo(
                click.style("Error: --card option requires a model ID", fg="red"),
                err=True,
            )
            return
        elif card and model:
            try:
                mc = ModelCard()
                model_metadata = mc.get(model, as_json=True)

                if not model_metadata:
                    click.echo(
                        click.style(
                            f"Error: No metadata found for model ID '{model}'", fg="red"
                        ),
                        err=True,
                    )
                    return
                click.echo(model_metadata)
            except Exception as e:
                click.echo(click.style(f"Error fetching model metadata: {e}", fg="red"))
            return
        else:
            mc = ModelCatalog(less=not more)

            if hub:
                catalog_table = mc.hub()
            else:
                catalog_table = mc.local()
                if not catalog_table.data:
                    click.echo(
                        click.style(
                            "No local models available. Please fetch a model by running 'ersilia fetch' command",
                            fg="red",
                        )
                    )
                    return
            if file_name is None:
                catalog = (
                    catalog_table.as_json() if as_json else catalog_table.as_table()
                )
            else:
                catalog_table.write(file_name)
                catalog = None

            click.echo(catalog)

    return catalog
