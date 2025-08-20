from ... import logger
from ...hub.delete.delete import ModelFullDeleter
from ..echo import echo

# from ...deletion import ModelFullDeleter  # adjust import as needed


def delete(model_id: str, verbose=False):
    """
    Deletes a specified model from local storage.

    Args:
        model_id (str): ID of the model to delete.

    Returns:
        str: Confirmation message on success or warning message on failure.

    Raises:
        RuntimeError: If the model cannot be deleted.
    """
    if verbose:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)

    md = ModelFullDeleter()
    can_delete, reason = md.can_be_deleted(model_id)
    # can_delete is bool, reason is message

    if can_delete:
        echo("Deleting model {0}".format(model_id))
        md.delete(model_id)
        echo(
            ":collision: Model {0} deleted successfully!".format(model_id),
            fg="green",
        )
    else:
        echo(
            f":person_tipping_hand: {reason}".format(model_id),
            fg="yellow",
        )
    return
