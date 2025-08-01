from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from ersilia.cli.commands.catalog import catalog_cmd


@pytest.fixture
def runner():
    return CliRunner()


@patch("ersilia.hub.content.catalog.ModelCatalog.airtable", return_value=None)
@patch("ersilia.hub.content.catalog.ModelCatalog.hub")
@patch("ersilia.hub.content.card.ModelCard.get", return_value={"model": "metadata"})
@pytest.mark.parametrize(
    "options, model",
    [
        (["--hub"], None),
    ],
)
def test_catalog_command(
    mock_model_card_get,
    mock_model_catalog_hub,
    mock_model_catalog_airtable,
    runner,
    options,
    model,
):
    mock_model_card_instance = MagicMock()

    if model:
        mock_model_card_instance.invoke.return_value.exit_code = 0
        mock_model_card_instance.invoke.return_value.output = (
            f"Fetching model {model}: \n👍 Model {model} executed successfully!\n"
        )

    result = runner.invoke(catalog_cmd(), options)
    assert (
        result.exit_code == 0
    ), f"Unexpected exit code: {result.exit_code}. Output: {result.output}"

    if "--card" in options and len(options) == 1:
        mock_model_card_get.assert_not_called()
    elif "--card" in options and model:
        mock_model_card_get.assert_called_once()
    elif "--browser" in options:
        mock_model_catalog_airtable.assert_called_once()
    elif "--hub" and "--local" in options:
        mock_model_catalog_hub.assert_not_called()
    elif "--hub" in options:
        mock_model_catalog_hub.assert_called_once()


if __name__ == "__main__":
    runner = CliRunner()
    test_catalog_command(None, None, runner, ["--local"], None)
