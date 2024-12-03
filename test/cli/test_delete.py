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
@pytest.mark.parametrize(
    "model",
    [
        (MODEL_ID)
    ]
)
def test_delete_model(
    mock_model_base,
    runner,
    model,
):
    mock_model_instance = MagicMock()
    mock_model_instance.model_id = model
    mock_model_base.return_value = mock_model_instance
    mock_model_instance.invoke.return_value.exit_code = 0
    mock_model_instance.invoke.return_value.output = (
        f"Deleting model {model}: \n👍 Model {model}\
          deleting cmd successfully executed!\n"
    )

    result = runner.invoke(delete_cmd(), [model])

    logger.info(result.output)

    assert (
        result.exit_code == 0
    ), f"Unexpected exit code: {result.exit_code}. Output: {result.output}"

if __name__ == "__main__":
    runner = CliRunner()
    test_delete_model(None, None, runner, MODEL_ID)
