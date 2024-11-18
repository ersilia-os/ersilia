import click
import csv

from ...hub.content.catalog import ModelCatalog


from . import ersilia_cli
from .. import echo
from ...hub.delete.delete import ModelFullDeleter
from ... import ModelBase


def delete_cmd():
    """Create delete command"""

    def _delete_model_by_id(model_id):
        md = ModelFullDeleter()
        if md.needs_delete(model_id):
            echo("Deleting model {0}".format(model_id))
            _delete(md, model_id)
            echo(
                ":collision: Model {0} deleted successfully!".format(model_id),
                fg="green",
            )
        else:
            echo(
                ":person_tipping_hand: Model {0} is not available locally. No delete is necessary".format(
                    model_id
                ),
                fg="yellow",
            )

    def _delete(md, model_id):
        md.delete(model_id)

    def _delete_all():
        """Function to delete all locally available models"""
        model_catalog = ModelCatalog()
        catalog_table = model_catalog.local()
        local_models = catalog_table.data if catalog_table else None
        if not local_models:
            echo(":person_tipping_hand: No models are available locally for deletion.", fg="yellow")
            return
        deleted_count = 0 
        for model_row in local_models:
            model_id = model_row[0]
            if _delete_model_by_id(model_id):
                deleted_count += 1
        echo(":thumbs_up: Completed the deletion of all locally available models!", fg="green")
        

    # Example usage:
    # 1. Delete a specific model: ersilia delete {MODEL}
    # 2. Delete all models: ersilia delete --all
    @ersilia_cli.command(
        short_help="Delete model from local computer",
        help="Delete model from local computer. The BentoML bundle is deleted, as well as the files stored in "
        "the EOS directory and the Pip-installed package",
    )
    @click.argument("model", required=False, type=click.STRING)
    @click.option("--all", is_flag=True, help="Delete all locally available models.")
    def delete(model, all):
        if all:
            _delete_all()
        elif model:
            model_id = ModelBase(model).model_id
            _delete_model_by_id(model_id)
        else:
            echo(":warning: Please specify a model to delete a model or use --all to delete all models.", fg="red")