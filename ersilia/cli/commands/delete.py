import rich_click as click

from ..echo import confirm, echo
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
        return md.delete(model_id)

    def _delete_model_by_id(model_id):
        from ...hub.delete.delete import ModelFullDeleter

        md = ModelFullDeleter()
        can_delete, reason = md.can_be_deleted(model_id)
        if can_delete:
            return _delete(md, model_id) is not False
        if "not available locally" in reason:
            echo(
                f"Model {model_id} is not available locally, so there is nothing to delete.",
                fg="yellow",
            )
        else:
            echo(reason, fg="red")
        return False

    def _delete_all():
        """Function to delete all locally available models"""
        from ...hub.content.catalog import ModelCatalog

        model_catalog = ModelCatalog()
        catalog_table = model_catalog.local()
        if not catalog_table:
            echo("No models are available locally.", fg="yellow")
            return
        local_models = catalog_table.data
        idx = catalog_table.columns.index("Identifier")
        if not local_models:
            echo("No models are available locally.", fg="yellow")
            return
        model_ids = [row[idx] for row in local_models]
        echo(
            "This will delete {0} model{1}: {2}.".format(
                len(model_ids), "" if len(model_ids) == 1 else "s", ", ".join(model_ids)
            ),
            fg="yellow",
        )
        if not confirm("Continue?", default=False):
            echo("Aborted. No models were deleted.")
            return
        deleted_count = 0
        for model_row in local_models:
            model_id = model_row[idx]
            try:
                if _delete_model_by_id(model_id):
                    deleted_count += 1
            except Exception as e:
                echo(f"Model {model_id} could not be deleted: {e}", fg="red")
        if deleted_count == len(local_models):
            echo(
                f"Deleted {deleted_count} model{'' if deleted_count == 1 else 's'}.",
                fg="green",
            )
        else:
            echo(
                f"Deleted {deleted_count} of {len(local_models)} models.",
                fg="red",
            )

    # Example usage:
    # 1. Delete a specific model: ersilia delete {MODEL}
    # 2. Delete all models: ersilia delete --all
    @ersilia_cli.command(
        short_help="Delete a model from this computer",
        help="Fully remove a model from the local computer. This includes the model files in the EOS directory, conda environment, Docker image and containers, pip package, and all associated database entries.",
    )
    @click.argument("model", required=False, type=click.STRING)
    @click.option("--all", is_flag=True, help="Delete all locally available models.")
    def delete(model, all):
        if all:
            _delete_all()
        elif model:
            from ... import ModelBase

            model_id = ModelBase(model).model_id
            _delete_model_by_id(model_id)
        else:
            raise click.UsageError(
                "Give a model to delete, or use --all to delete every local model."
            )

    return delete
