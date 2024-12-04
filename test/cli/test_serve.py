import pytest
from unittest import TestCase
from unittest.mock import patch, MagicMock, AsyncMock
from click.testing import CliRunner
from ersilia.core.model import ErsiliaModel
from ersilia.cli.commands.serve import serve_cmd
from ersilia.store.utils import store_has_model

URL = "http://localhost"
MODEL_ID = "eos3b5e"


@pytest.fixture
def mock_set_apis():
    with patch.object(
        ErsiliaModel, 
        "_set_apis", 
        return_value=None
        ) as mock_set_apis:
        yield mock_set_apis

@pytest.fixture
def mock_serve():
    with patch.object(
        ErsiliaModel, 
        "serve", 
        return_value=None
        ) as mock_serve_:
        yield mock_serve_



@patch("ersilia.core.model.ErsiliaModel")
@patch("ersilia.store.utils.store_has_model", return_value=False)
def test_serve_cmd(
    mock_store_has_model, 
    mock_ersilia_model,
    mock_set_apis,
    mock_serve
):
    runner = CliRunner()
    mock_mdl_instance = MagicMock()
    mock_mdl_instance.is_valid.return_value = True
    mock_mdl_instance.url = URL
    mock_mdl_instance.model_id = MODEL_ID
    mock_mdl_instance.slug = "molecular-weight"
    mock_mdl_instance.pid = 1234
    mock_mdl_instance.scl = "docker"
    mock_mdl_instance.output_source = "LOCAL_ONLY"
    mock_mdl_instance.get_apis.return_value = ["run"]

    mock_ersilia_model.return_value = mock_mdl_instance

    result = runner.invoke(serve_cmd(), [MODEL_ID, "--docker"])

    assert result.exit_code == 0
    assert mock_serve.called