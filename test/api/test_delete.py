import pytest
from unittest.mock import MagicMock, patch

# from ersilia.api import ErsiliaAPI
import ersilia.api.commands.delete as delete_mod

MODEL_ID = "eos3b5e"
REASON = "model is in use and cannot be deleted"

@patch.object(delete_mod, "logger")
@patch.object(delete_mod, "echo")
@patch.object(delete_mod, "ModelFullDeleter")
def test_delete_success(mock_deleter_cls, mock_echo, mock_logger):
    # Arrange: can_be_deleted returns True
    mock_md = MagicMock()
    mock_md.can_be_deleted.return_value = (True, "")
    mock_deleter_cls.return_value = mock_md

    # Act: call with verbose=False (default)
    delete_mod.delete(MODEL_ID)

    # Assert: verbosity was set to 0
    mock_logger.set_verbosity.assert_called_once_with(0)
    # Assert: can_be_deleted was called with our model_id
    mock_md.can_be_deleted.assert_called_once_with(MODEL_ID)
    # Assert: deletion path: echo “Deleting…”, then md.delete(), then success echo
    mock_echo.assert_any_call(f"Deleting model {MODEL_ID}")
    mock_md.delete.assert_called_once_with(MODEL_ID)
    mock_echo.assert_any_call(
        f":collision: Model {MODEL_ID} deleted successfully!", fg="green"
    )

@patch.object(delete_mod, "logger")
@patch.object(delete_mod, "echo")
@patch.object(delete_mod, "ModelFullDeleter")
def test_delete_failure(mock_deleter_cls, mock_echo, mock_logger):
    # Arrange: can_be_deleted returns False with a reason
    mock_md = MagicMock()
    mock_md.can_be_deleted.return_value = (False, REASON)
    mock_deleter_cls.return_value = mock_md

    # Act: call with verbose=True
    delete_mod.delete(MODEL_ID, verbose=True)

    # Assert: verbosity was set to 1
    mock_logger.set_verbosity.assert_called_once_with(1)
    # Assert: can_be_deleted was called correctly
    mock_md.can_be_deleted.assert_called_once_with(MODEL_ID)
    # Assert: delete() was NOT called
    mock_md.delete.assert_not_called()
    # Assert: echo was called once with the reason
    mock_echo.assert_called_once_with(f":person_tipping_hand: {REASON}", fg="yellow")
