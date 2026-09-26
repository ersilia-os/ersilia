import sys

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
            return None
        echo(reason, fg="red")
        return False

    def _close_if_served(model_ids):
        # A model served in this terminal is closed before it is deleted, so
        # its container and session record do not outlive it.
        from ... import ErsiliaModel
        from ...core.session import Session
        from ...utils.session import deregister_model_session

        session = Session(config_json=None)
        served = _served_model()
        if served not in model_ids:
            return
        mdl = ErsiliaModel(served, service_class=session.current_service_class())
        mdl.close()
        deregister_model_session(served)
        echo(f"Model {served} closed.", fg="green")

    def _served_model():
        from ...core.session import Session

        try:
            return Session(config_json=None).current_model_id()
        except Exception:
            # A half-written session record: nothing usable is served.
            return None

    def _served_here(model_id):
        return _served_model() == model_id

    def _delete_all():
        """Function to delete all locally available models"""
        from ...hub.content.catalog import ModelCatalog

        model_catalog = ModelCatalog()
        catalog_table = model_catalog.local()
        model_ids = []
        if catalog_table and catalog_table.data:
            idx = catalog_table.columns.index("Identifier")
            model_ids = [row[idx] for row in catalog_table.data]
        # Remains of models not fully fetched or deleted go too.
        model_ids += [m for m in model_catalog.leftovers if m not in model_ids]
        if not model_ids:
            echo("No models are available locally.", fg="yellow")
            return
        echo(
            "This will delete {0} model{1}: {2}.".format(
                len(model_ids), "" if len(model_ids) == 1 else "s", ", ".join(model_ids)
            ),
            fg="yellow",
        )
        if not confirm("Continue?", default=False):
            echo("Aborted. No models were deleted.")
            return
        _close_if_served(model_ids)
        deleted_count = 0
        for model_id in model_ids:
            try:
                if _delete_model_by_id(model_id):
                    deleted_count += 1
            except Exception as e:
                echo(f"Model {model_id} could not be deleted: {e}", fg="red")
        if deleted_count == len(model_ids):
            echo(
                f"Deleted {deleted_count} model{'' if deleted_count == 1 else 's'}.",
                fg="green",
            )
        else:
            echo(
                f"Deleted {deleted_count} of {len(model_ids)} models.",
                fg="red",
            )
            sys.exit(1)

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
            if _served_here(model_id):
                if not confirm(
                    f"Model {model_id} is being served in this terminal. Close and delete it?",
                    default=False,
                ):
                    echo(f"Aborted. Model {model_id} was not deleted.")
                    return
                _close_if_served([model_id])
            if _delete_model_by_id(model_id) is False:
                sys.exit(1)
        else:
            raise click.UsageError(
                "Give a model to delete, or use --all to delete every local model."
            )

    return delete
