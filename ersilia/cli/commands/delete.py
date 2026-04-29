import rich_click as click

from ... import ModelBase
from ...hub.content.catalog import ModelCatalog
from ...hub.delete.delete import ModelFullDeleter
from .. import echo
from . import ersilia_cli


def delete_cmd():
    """
    Deletes a specified model.

    This command allows users to delete a specified model from the local storage.

    Returns
    -------
    function
        The delete command function to be used by the CLI and for testing in the pytest.

    Examples
    --------
    .. code-block:: console

        Delete a specific model:
        $ ersilia delete <model_id>

        Delete all models:
        $ ersilia delete --all
    """

    def _delete(md, model_id):
        md.delete(model_id)

    def _delete_model_by_id(model_id):
        md = ModelFullDeleter()
        can_delete, reason = md.can_be_deleted(model_id)
        if can_delete:
            echo("Deleting model {0}".format(model_id), harmonize=False)
            _delete(md, model_id)
            return True
        else:
            if "not available locally" in reason:
                echo(f"Model {model_id} is not available locally.", fg="red")
            else:
                echo(f"{reason}", fg="red")
            return False

    def _delete_all():
        """Function to delete all locally available models"""
        model_catalog = ModelCatalog()
        catalog_table = model_catalog.local()
        local_models = catalog_table.data if catalog_table else None
        idx = catalog_table.columns.index("Identifier")
        if not local_models:
            echo(
                "No local models available.",
                fg="red",
            )
            return
        model_ids = [row[idx] for row in local_models]
        echo(
            "This will delete the following models: {0}".format(", ".join(model_ids)),
            fg="yellow",
        )
        if not click.confirm("Are you sure you want to proceed?", default=False):
            echo("Aborted.", fg='green')
            return
        deleted_count = 0
        for model_row in local_models:
            model_id = model_row[idx]
            try:
                if _delete_model_by_id(model_id):
                    deleted_count += 1
            except Exception as e:
                echo(
                    f":warning: Error deleting model {model_id}: {e}",
                    fg="red",
                )
        if deleted_count is len(local_models):
            echo(
                f":thumbs_up: Completed the deletion of {deleted_count} locally available models!",
                fg="green",
            )
        else:
            echo(
                f":thumbs_down: Failed to delete all models. Deleted {deleted_count} out of {len(local_models)} models.",
                fg="red",
            )

    # Example usage:
    # 1. Delete a specific model: ersilia delete {MODEL}
    # 2. Delete all models: ersilia delete --all
    @ersilia_cli.command(
        short_help="Delete model from local computer",
        help="Fully remove a model from the local computer. This includes the model files in the EOS directory, conda environment, Docker image and containers, pip package, and all associated database entries.",
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
            echo(
                ":warning: Please specify a model to delete a model or use --all to delete all models.",
                fg="red",
            )

    return delete
