from unittest.mock import MagicMock, patch
import pytest

from ersilia.api.commands.serve import serve
from ersilia.core.model import ErsiliaModel
from ersilia.hub.fetch.register.standard_example import ModelStandardExample #not on chat
from ersilia.store.utils import OutputSource #from chat

URL = "http://localhost"
MODEL_ID = "eos3b5e"
SESSION_DIR = "/mock/session"
SCL = "pulled_docker"

# Fixture to mock the ErsiliaModel used inside the serve() function
@pytest.fixture
def mock_ersilia_model():
    # Patch ErsiliaModel inside the API's serve module
    with patch("ersilia.api.commands.serve.ErsiliaModel") as MockModel:
        # Create a mock model instance with the expected attributes/methods
        mock_instance = MagicMock()
        mock_instance.is_valid.return_value = True
        mock_instance.url = URL  # URL that would be returned by the running model
        mock_instance.session._session_dir = SESSION_DIR  # Fake session directory
        mock_instance.scl = SCL  # Service class label (e.g., pulled_docker)
        mock_instance.get_apis.return_value = ["run"]  # API exposes only 'run'

        # Set the patched class to return the mock instance
        MockModel.return_value = mock_instance
        yield mock_instance

# Fixture to mock Redis cache setup logic (SetupRedis class)
@pytest.fixture
def mock_redis_setup():
    with patch("ersilia.api.commands.serve.SetupRedis") as MockRedis:
        yield MockRedis

# Fixture to mock session registration (register_model_session)
@pytest.fixture
def mock_register_model_session():
    with patch("ersilia.api.commands.serve.register_model_session") as mock_reg:
        yield mock_reg

# The actual test case for the API's serve() function
def test_serve_api_basic(
    mock_ersilia_model, mock_redis_setup, mock_register_model_session
):
    # Call the API version of serve with standard arguments
    url, session_dir, scl = serve_api(
        model=MODEL_ID,
        port=None,
        track=False,
        tracking_use_case="local",
        enable_local_cache=True,
        local_cache_only=False,
        cloud_cache_only=False,
        cache_only=False,
        max_cache_memory_frac=0.3,
        verbose=False,
    )

    # ========== ASSERTIONS ==========

    # Check that the return values match the mocked model attributes
    assert url == URL
    assert session_dir == SESSION_DIR
    assert scl == SCL

    # Ensure that serve() was called on the mock model
    mock_ersilia_model.serve.assert_called_once()

    # Ensure register_model_session was called with the correct args
    mock_register_model_session.assert_called_once_with(MODEL_ID, SESSION_DIR)

    # Ensure SetupRedis was initialized with expected values
    mock_redis_setup.assert_called_once_with(True, 0.3)