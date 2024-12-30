import pytest
from unittest import TestCase
from unittest.mock import patch, MagicMock, AsyncMock
from click.testing import CliRunner
from ersilia.core.model import ErsiliaModel
from ersilia.cli.commands.close import close_cmd
from ersilia.core.session import Session

MODEL_ID = "eos3b5e"


@pytest.fixture
def mock_close():
    with patch.object(ErsiliaModel, "close", return_value=None) as mock_close_:
        yield mock_close_


@pytest.fixture
def mock_session():
    with (
        patch.object(Session, "current_model_id", return_value=MODEL_ID),
        patch.object(Session, "current_service_class", return_value="pulled_docker"),
        patch.object(Session, "tracking_status", return_value=False),
        patch.object(Session, "current_output_source", return_value="LOCAL_ONLY"),
    ):
        yield


def test_close_cmd(mock_close, mock_session):
    runner = CliRunner()

    result = runner.invoke(close_cmd())

    assert result.exit_code == 0
    # assert mock_fetcher.called
    assert mock_close.called
