# Adjust the import paths as needed based on your project structure
from ersilia.hub.content.catalog import ModelCatalog
from ersilia.hub.delete.delete import ModelFullDeleter


def _delete_model_by_id(model_id):
    deleter = ModelFullDeleter(model_id)
    deleter.delete()


def delete_all_models() -> dict:
    """
    Deletes all locally available models.

    Returns:
        dict: A summary of the deletion process, including total models deleted
              and any errors encountered.
    """
    model_catalog = ModelCatalog()
    catalog_table = model_catalog.local()
    local_models = catalog_table.data if catalog_table else None

    if not local_models:
        print(" No models are available locally for deletion.")
        return {
            "deleted_count": 0,
            "errors": [],  # or None?
        }

    idx = catalog_table.columns.index("Identifier")  # diff placement
    deleted_count = 0
    errors = []

    for model_row in local_models:
        model_id = model_row[idx]
        try:
            _delete_model_by_id(model_id)
            deleted_count += 1
        except Exception as e:
            error_msg = f"Error deleting model {model_id}: {e}"
            print(error_msg)
            errors.append(error_msg)

    print(f" Completed the deletion of {deleted_count} locally available models!")

    return {"deleted_count": deleted_count, "errors": errors}
