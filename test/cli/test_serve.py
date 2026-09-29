from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from ersilia.cli.commands.serve import serve_cmd
from ersilia.core.model import ErsiliaModel
from ersilia.hub.fetch.register.standard_example import ModelStandardExample

URL = "http://localhost"
MODEL_ID = "eos3b5e"


@pytest.fixture
def mock_std_example():
    with patch.object(ModelStandardExample, "run", return_value=None) as mock_run:
        yield mock_run


@pytest.fixture
def mock_serve():
    with patch.object(ErsiliaModel, "serve", return_value=None) as mock_serve_:
        yield mock_serve_


@patch("ersilia.core.session.Session.get", return_value={})
@patch("ersilia.utils.session.register_model_session")
@patch("ersilia.core.model.ErsiliaModel")
@patch("ersilia.store.utils.store_has_model", return_value=False)
def test_serve_cmd(
    mock_store_has_model,
    mock_ersilia_model,
    mock_register,
    mock_session_get,
    mock_serve,
    mock_std_example,
):
    # register_model_session is mocked so the test never writes to the real
    # ~/eos/models.json (a mocked model id there corrupted the file).
    """Verify that the serve command exits successfully and calls ErsiliaModel.serve for the given model ID."""
    runner = CliRunner()
    mock_mdl_instance = MagicMock()
    mock_mdl_instance.is_valid.return_value = True
    mock_mdl_instance.url = URL
    mock_mdl_instance.model = MODEL_ID
    mock_mdl_instance.service_class = MODEL_ID
    mock_mdl_instance.slug = "molecular-weight"
    mock_mdl_instance.pid = 1234
    mock_mdl_instance.scl = "pulled_docker"
    mock_mdl_instance.output_source = "LOCAL_ONLY"
    mock_mdl_instance.get_apis.return_value = ["run"]
    mock_mdl_instance.model_id = MODEL_ID
    mock_mdl_instance.info.return_value = {"docker_tag": "latest", "card": {}}

    mock_ersilia_model.return_value = mock_mdl_instance

    # The command imports ErsiliaModel from the ersilia package, which may
    # already hold the real class; patch it there too.
    with patch("ersilia.ErsiliaModel", mock_ersilia_model, create=True):
        result = runner.invoke(serve_cmd(), [MODEL_ID])

    assert result.exit_code == 0
    # serve is called on the mocked instance or the patched method, depending
    # on which ErsiliaModel the command imports first.
    assert mock_serve.called or mock_mdl_instance.serve.called
