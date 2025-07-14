from unittest.mock import MagicMock, patch
import pytest

from ersilia.api.create_api import ErsiliaAPI
from ersilia.hub.delete.delete import ModelFullDeleter

MODEL = "eos3b5e"
DUMMY_MODEL = ["eosxxxx"]


@pytest.fixture
def can_be_deleted():
    with patch.object(
        ModelFullDeleter, "can_be_deleted", return_value=(True, "")
    ) as can_be_deleted:
        yield can_be_deleted


@patch("ersilia.ModelBase")
@patch("ersilia.hub.delete.delete.ModelFullDeleter")
@patch("ersilia.cli.echo")
def test_delete_specific_model(mock_deleter, mock_modelbase, can_be_deleted):
    mock_modelbase_instance = MagicMock()
    mock_modelbase_instance.model_id = MODEL  #
    mock_modelbase.return_value = mock_modelbase_instance

    mock_deleter_instance = MagicMock()
    mock_deleter_instance.can_be_deleted.return_value = (True, "")
    mock_deleter_instance.delete.return_value = None
    mock_deleter.return_value = mock_deleter_instance

    api = ErsiliaAPI(MODEL)
    result = api.delete()

    mock_deleter_instance.delete.assert_called_once()
    assert result is None or result == "Deleted", "API did not behave as expected."


@patch("ersilia.hub.content.catalog.ModelCatalog")
@patch("ersilia.hub.delete.delete.ModelFullDeleter")
@patch("ersilia.cli.echo")
def test_delete_all_models(mock_deleter, mock_catalog):
    mock_catalog_instance = MagicMock()
    mock_catalog_instance.local.return_value.data = [[MODEL], [DUMMY_MODEL]]
    mock_catalog_instance.local.return_value.columns = ["Identifier"]
    mock_catalog.return_value = mock_catalog_instance

    mock_deleter_instance = MagicMock()
    mock_deleter_instance.can_be_deleted.return_value = (True, "")
    mock_deleter_instance.delete.return_value = None
    mock_deleter.return_value = mock_deleter_instance

    for model_id in [MODEL, DUMMY_MODEL]:
        api = ErsiliaAPI(model_id)
        api.delete()

    assert mock_deleter_instance.delete.call_count == 2


@patch("ersilia.hub.content.catalog.ModelCatalog")
@patch("ersilia.cli.echo")
def test_no_models_available(mock_echo, mock_catalog):
    runner = CliRunner()

    mock_catalog_instance = MagicMock()
    mock_catalog_instance.local.return_value = None
    mock_catalog.return_value = mock_catalog_instance

    result = runner.invoke(delete(), ["delete", "--all"])

    assert (
        result.exit_code == 0
    ), f"Unexpected exit code: {result.exit_code}. Output: {result.output}"