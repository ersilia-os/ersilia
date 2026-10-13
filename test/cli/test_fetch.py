from unittest.mock import MagicMock, patch

import pytest
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
        (
            "molecular-weight",
            MODEL_ID,
            ["--from_apptainer"],
        ),  # Test with --from_apptainer flag
        (
            "molecular-weight",
            MODEL_ID,
            ["--from_apptainer", "--version", "v1"],
        ),  # --version applies to Apptainer too
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
    """Verify that fetching a known model exits successfully, with and without the --from_dockerhub flag."""
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
    """Verify that fetching an unrecognised model ID exits with code 1."""
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


@patch("ersilia.core.modelbase.ModelBase")
@patch(
    "ersilia.hub.fetch.fetch.ModelFetcher.__init__",
    return_value=None,
)
@patch(
    "ersilia.hub.fetch.fetch.ModelFetcher.fetch",
    return_value=FetchResult(True, "Model fetched successfully"),
)
def test_fetch_from_apptainer_is_passed_on(
    mock_fetch, mock_init, mock_model_base, runner
):
    """Verify that --from_apptainer reaches the fetcher, and DockerHub is not forced."""
    mock_model_base.return_value = MagicMock(model_id=MODEL_ID, slug="slug")
    result = runner.invoke(fetch_cmd(), [MODEL_ID, "--from_apptainer"])
    assert result.exit_code == 0, result.output
    kwargs = mock_init.call_args.kwargs
    assert kwargs["force_from_apptainer"] is True
    assert kwargs["force_from_dockerhub"] is False


def test_fetch_from_apptainer_and_dockerhub_conflict(runner):
    """Verify that asking for two sources is refused."""
    result = runner.invoke(
        fetch_cmd(), [MODEL_ID, "--from_apptainer", "--from_dockerhub"]
    )
    assert result.exit_code != 0
    assert "Choose only one source" in result.output
