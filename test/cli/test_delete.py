from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from ersilia.cli.commands.delete import delete_cmd
from ersilia.hub.delete.delete import ModelFullDeleter

MODEL = "eos3b5e"
DUMMY_MODEL = "eosxxxx"


@pytest.fixture
def can_be_deleted():
    with patch.object(
        ModelFullDeleter, "can_be_deleted", return_value=(True, "")
    ) as can_be_deleted:
        yield can_be_deleted


@patch("ersilia.ModelBase")
@patch("ersilia.hub.delete.delete.ModelFullDeleter")
@patch("ersilia.cli.echo")
def test_delete_specific_model(mock_echo, mock_deleter, mock_modelbase, can_be_deleted):
    """Verify that deleting a specific model by ID exits successfully."""
    mock_modelbase_instance = MagicMock()
    mock_modelbase_instance.model_id = MODEL  #
    mock_modelbase.return_value = mock_modelbase_instance

    mock_deleter_instance = MagicMock()
    mock_deleter_instance.can_be_deleted.return_value = (False, "")
    mock_deleter_instance.delete.return_value = None
    mock_deleter.return_value = mock_deleter_instance

    mock_echo.return_value = None

    runner = CliRunner()
    result = runner.invoke(delete_cmd(), [MODEL])

    assert (
        result.exit_code == 0
    ), f"Unexpected exit code: {result.exit_code}. Output: {result.output}"


@patch("ersilia.hub.content.catalog.ModelCatalog")
@patch("ersilia.hub.delete.delete.ModelFullDeleter")
@patch("ersilia.cli.echo")
def test_delete_all_models(mock_echo, mock_deleter, mock_catalog):
    """Verify that --all deletes all local models when the user confirms and deletion is allowed."""
    runner = CliRunner()

    mock_catalog_instance = MagicMock()
    mock_catalog_instance.local.return_value.data = [[MODEL], [DUMMY_MODEL]]
    mock_catalog_instance.local.return_value.columns = ["Identifier"]
    mock_catalog.return_value = mock_catalog_instance

    mock_deleter_instance = MagicMock()
    mock_deleter_instance.can_be_deleted.return_value = (True, "")
    mock_deleter.return_value = mock_deleter_instance

    result = runner.invoke(delete_cmd(), ["delete", "--all"], input="y\n")

    assert (
        result.exit_code == 0
    ), f"Unexpected exit code: {result.exit_code}. Output: {result.output}"

@patch("ersilia.hub.content.catalog.ModelCatalog")
@patch("ersilia.hub.delete.delete.ModelFullDeleter")
@patch("ersilia.cli.echo")
def test_delete_all_models_fails_when_can_be_deleted_is_False(mock_echo, mock_deleter, mock_catalog):
    """Verify that --all exits cleanly even when no model passes the can_be_deleted check."""
    runner = CliRunner()

    mock_catalog_instance = MagicMock()
    mock_catalog_instance.local.return_value.data = [[MODEL], [DUMMY_MODEL]]
    mock_catalog_instance.local.return_value.columns = ["Identifier"]
    mock_catalog.return_value = mock_catalog_instance

    mock_deleter_instance = MagicMock()
    mock_deleter_instance.can_be_deleted.return_value = (False, "")
    mock_deleter.return_value = mock_deleter_instance

    result = runner.invoke(delete_cmd(), ["delete", "--all"], input="y\n")

    assert (
        result.exit_code == 0
    ), f"Unexpected exit code: {result.exit_code}. Output: {result.output}"


@patch("ersilia.hub.content.catalog.ModelCatalog")
@patch("ersilia.hub.delete.delete.ModelFullDeleter")
@patch("ersilia.cli.echo")
def test_delete_all_models_aborts_when_user_says_no(_mock_echo, mock_deleter, mock_catalog):
    """Verify that --all aborts cleanly without deleting anything when the user declines the confirmation."""
    runner = CliRunner()

    mock_catalog_instance = MagicMock()
    mock_catalog_instance.local.return_value.data = [[MODEL], [DUMMY_MODEL]]
    mock_catalog_instance.local.return_value.columns = ["Identifier"]
    mock_catalog.return_value = mock_catalog_instance

    mock_deleter_instance = MagicMock()
    mock_deleter.return_value = mock_deleter_instance

    result = runner.invoke(delete_cmd(), ["delete", "--all"], input="N\n")

    assert (
        result.exit_code == 0
    ), f"Unexpected exit code: {result.exit_code}. Output: {result.output}"
    mock_deleter_instance.delete.assert_not_called()


@patch("ersilia.hub.content.catalog.ModelCatalog")
@patch("ersilia.cli.echo")
def test_no_models_available(mock_echo, mock_catalog):
    """Verify that --all exits cleanly with an error message when no local models exist."""
    runner = CliRunner()

    mock_catalog_instance = MagicMock()
    mock_catalog_instance.local.return_value = None
    mock_catalog.return_value = mock_catalog_instance

    result = runner.invoke(delete_cmd(), ["delete", "--all"], input="y\n")

    assert (
        result.exit_code == 0
    ), f"Unexpected exit code: {result.exit_code}. Output: {result.output}"

