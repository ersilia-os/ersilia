import pytest
from unittest.mock import patch, MagicMock
from click.testing import CliRunner
from ersilia.cli.commands.delete import delete_cmd
from ersilia.utils.logging import logger

@pytest.fixture
def runner():
    return CliRunner()

MODEL_ID = "eos3b5e"
@patch("ersilia.core.modelbase.ModelBase")
@patch("ersilia.hub.delete.delete.ModelFullDeleter.needs_delete", return_value=None)
@pytest.mark.parametrize(
    "model",
    [
        (MODEL_ID)
    ]
)
def test_delete_model(
    mock_delete,
    mock_model_base,
    runner,
    model,
):
    mock_model_instance = MagicMock()
    mock_model_instance.model_id = model
    mock_model_base.return_value = mock_model_instance
    mock_model_instance.invoke.return_value.exit_code = 0
    mock_model_instance.invoke.return_value.output = (
        f"Fetching model {model}: \n👍 Model {model} deleted successfully!\n"
    )

    result = runner.invoke(delete_cmd(), [model])

    logger.info(result.output)

    assert (
        result.exit_code == 0
    ), f"Unexpected exit code: {result.exit_code}. Output: {result.output}"

    mock_delete.assert_called_once()

if __name__ == "__main__":
    runner = CliRunner()
    # Directly execute the test without pytest
    test_delete_model(None, None, runner, MODEL_ID)
