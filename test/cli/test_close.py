import pytest

from click.testing import CliRunner

from ersilia.cli.commands.close import close_cmd
from ersilia.core.model import ErsiliaModel
from ersilia.core.session import Session
from unittest.mock import patch, AsyncMock

MODEL_ID = "eos3b5e"
URL = "http://localhost"
PORT = 8001
API_NAME = "run"


@pytest.fixture
def mock_close():
    with patch.object(ErsiliaModel, "close", return_value=None) as mock_close_:
        yield mock_close_


@pytest.fixture
def mock_fetcher():
    with patch(
        "ersilia.hub.fetch.fetch.ModelFetcher.fetch", new_callable=AsyncMock
    ) as mock_fetch:
        yield mock_fetch


@pytest.fixture
def mock_set_apis():
    with patch.object(ErsiliaModel, "_set_apis", return_value=None) as mock_set_apis:
        yield mock_set_apis


@pytest.fixture
def mock_convn_api_get_apis():
    def mock_get_api_side_effect():
        return [API_NAME]

    with patch.object(
        ErsiliaModel, "get_apis", side_effect=mock_get_api_side_effect
    ) as mock_get_apis:
        yield mock_get_apis


@pytest.fixture
def mock_get_url():
    with patch.object(
        ErsiliaModel, "_get_url", return_value=f"{URL}:{PORT}"
    ) as mock_url:
        yield mock_url


@pytest.fixture
def mock_session():
    with (
        patch.object(Session, "current_model_id", return_value=MODEL_ID),
        patch.object(Session, "current_service_class", return_value="pulled_docker"),
        patch.object(Session, "tracking_status", return_value=False),
        patch.object(Session, "current_output_source", return_value="LOCAL_ONLY"),
    ):
        yield


def test_close_cmd(
    mock_close,
    mock_fetcher,
    mock_session,
    mock_set_apis,
    mock_convn_api_get_apis,
    mock_get_url,
):
    runner = CliRunner()

    result = runner.invoke(close_cmd())

    assert result.exit_code == 0
    assert mock_close.called
