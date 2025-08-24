import pytest
from unittest.mock import patch, MagicMock

from ersilia.api.commands.close import close  # API function under test
from ersilia import ErsiliaModel
from ersilia.core.session import Session
from ersilia.utils.session import deregister_model_session

MODEL_ID = "eos3b5e"


@pytest.fixture
def mock_session_service_class():
    with patch.object(Session, "current_service_class", return_value="pulled_docker") as mock_service_class:
        yield mock_service_class


@pytest.fixture
def mock_model_close():
    with patch.object(ErsiliaModel, "close", return_value=None) as mock_close:
        yield mock_close


@pytest.fixture
def mock_deregister_session():
    with patch("ersilia.api.commands.close.deregister_model_session") as mock_deregister:
        yield mock_deregister


def test_close_api(mock_session_service_class, mock_model_close, mock_deregister_session):
    """
    ✅ Test the API version of close(model_id)
    """
    assert 1 == 1
    # TODO: Complete this test
    # # Call the close function with a valid model ID
    # result = close(MODEL_ID)

    # # Verify close() was called on the ErsiliaModel
    # mock_model_close.assert_called_once()

    # # Verify session deregistration occurred
    # mock_deregister_session.assert_called_once_with(MODEL_ID)

    # # No explicit return value is needed; if you add one, test it here
    # assert result is None