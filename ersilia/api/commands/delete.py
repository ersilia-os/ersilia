from ... import ModelBase
from ...hub.content.catalog import ModelCatalog
from ...hub.delete.delete import ModelFullDeleter

#from ...deletion import ModelFullDeleter  # adjust import as needed

def delete(model_id: str) -> str:
    """
    Deletes a specified model from local storage.

    Args:
        model_id (str): ID of the model to delete.

    Returns:
        str: Confirmation message on success or warning message on failure.

    Raises:
        RuntimeError: If the model cannot be deleted.
    """
    md = ModelFullDeleter()
    can_delete, reason = md.can_be_deleted(model_id)
    #can_delete is bool, reason is message

    if can_delete:
        print(f"Deleting model {model_id}...")
        md.delete(model_id)
        return f" Model {model_id} deleted successfully!"
    else:
        raise RuntimeError(f" {reason}")
