import pytest
from unittest.mock import patch, MagicMock
from click.testing import CliRunner
from ersilia.cli.commands.fetch import fetch_cmd
from ersilia.hub.fetch.fetch import FetchResult
from ersilia.utils.logging import logger

MODEL_ID = "eos3b5e"


@pytest.fixture
def runner():
    return CliRunner()


@patch("ersilia.core.modelbase.ModelBase")
@patch(
    "ersilia.hub.fetch.fetch.ModelFetcher.fetch",
    return_value=FetchResult(True, "Model fetched successfully."),
)
@pytest.mark.parametrize(
    "slug, model, flags",
    [
        ("molecular-weight", MODEL_ID, []),  # Test with no flags
        (
            "molecular-weight",
            MODEL_ID,
            ["--from_dockerhub"],
        ),  # Test with --from_dockerhub flag
    ],
)
def test_fetch_multiple_model(
    mock_fetch,
    mock_model_base,
    runner,
    slug,
    model,
    flags,
):
    """
    Fetch multiple model based on the mocking
    ============================================
    This codes creates a mock instance of a model for testing purposes.
    It sets the model's ID and slug attributes to simulate a real model,
    and configures the ModelBase mock to return this instance.
    Additionally, it mocks the behavior of the model's invoke method to
    simulate a successful command execution, specifying an exit code of 0
    and providing a formatted output string that mimics the expected
    success message when fetching the model. This allows for testing
    without the need for actual model fetching, ensuring that tests can
    run independently and reliably.
    """
    if flags is None:
        flags = []

    mock_model_instance = MagicMock()
    mock_model_instance.model_id = model
    mock_model_instance.slug = slug
    mock_model_base.return_value = mock_model_instance
    mock_model_instance.invoke.return_value.exit_code = 0
    mock_model_instance.invoke.return_value.output = (
        f"Fetching model {model}: {slug}\n👍 Model {model} fetched successfully!\n"
    )

    result = runner.invoke(fetch_cmd(), [model] + flags)

    logger.info(result.output)

    # Ensures that the click command run and exited corectly
    assert (
        result.exit_code == 0
    ), f"Unexpected exit code: {result.exit_code}. Output: {result.output}"

    # Ensure the fetch method was called with expected arguments
    # This will actually ensures the mocking did actually correctly implemented
    mock_fetch.assert_called_once()


@patch("ersilia.core.modelbase.ModelBase")
@patch("ersilia.hub.fetch.fetch.ModelFetcher.fetch", return_value=None)
def test_fetch_unknown_model(
    mock_fetch,
    mock_model_base,
    runner,
    slug="random-slug",
    model="xeos3111",  # something different
    flags=["--from_dockerhub"],
):
    if flags is None:
        flags = []

    mock_model_instance = MagicMock()
    mock_model_instance.model_id = model
    mock_model_instance.slug = slug
    mock_model_base.return_value = mock_model_instance
    mock_model_instance.invoke.return_value.exit_code = 1
    mock_model_instance.invoke.return_value.output = (
        f"Fetching model {model}: {slug}\nModel not found!\n"
    )

    result = runner.invoke(fetch_cmd(), [model] + flags)
    # This will create mess in the terminal (Something went wrong with Ersilia)
    # Which is expected as well with this exeption: Ersilia exception class: InvalidModelIdentifierError

    logger.info(result.output)

    # Ensures that the click command run and exited corectly
    assert (
        result.exit_code == 1
    ), f"Unexpected exit code: {result.exit_code}. Output: {result.output}"


if __name__ == "__main__":
    runner = CliRunner()
    # Directly execute the test without pytest
    test_fetch_multiple_model(None, None, runner, "molecular-weight", MODEL_ID, [])
    test_fetch_multiple_model(
        None, None, runner, "molecular-weight", MODEL_ID, ["--from_dockerhub"]
    )
    test_fetch_multiple_model(
        None, None, runner, "molecular-weight", MODEL_ID, ["--from_github"]
    )
    test_fetch_unknown_model(None, None, runner)
