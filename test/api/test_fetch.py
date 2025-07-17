from unittest.mock import MagicMock, patch
import pytest

from ersilia.api.commands.fetch import fetch
from ersilia.hub.fetch.fetch import FetchResult
from ersilia.utils.logging import logger

MODEL_ID = "eos3b5e"
SLUG = "molecular-weight" #chat addition

# Patch the asynchronous fetch inside the ModelFetcher
@pytest.fixture
def mock_model_fetcher():
    with patch("ersilia.api.commands.fetch.ModelFetcher.fetch") as mock_fetch:
        mock_fetch.return_value = FetchResult(True, "Model fetched successfully.")
        yield mock_fetch

# Patch ModelBase to return a mocked model with model_id and slug
@pytest.fixture
def mock_model_base():
    with patch("ersilia.api.commands.fetch.ModelBase") as MockModelBase:
        mock_instance = MagicMock()
        mock_instance.model_id = MODEL_ID
        mock_instance.slug = SLUG
        MockModelBase.return_value = mock_instance
        yield mock_instance


def test_fetch_api_success(mock_model_fetcher, mock_model_base):
    """
    ✅ Test successful fetch using the API function.
    """

    result = fetch(
        model=MODEL_ID,
        overwrite=True,
        from_dir=None,
        from_github=False,
        from_dockerhub=True,
        version="latest",
        from_s3=False,
        from_hosted=False,
        hosted_url=None,
        verbose=True,
    )

    # Log output for visibility during test runs
    logger.info(f"✅ Fetch succeeded: {mock_model_fetcher.return_value.reason}")

    # Assert result is success
    assert result is not False

    # Validate the mock fetch call
    mock_model_fetcher.assert_called_once()
    mock_model_base.assert_called_once()


def test_fetch_api_failure():
    """
    ❌ Test failed fetch using the API function.
    """

    with patch("ersilia.api.commands.fetch.ModelBase") as mock_model_base, \
         patch("ersilia.api.commands.fetch.ModelFetcher.fetch") as mock_fetch:

        # Simulate an invalid model
        model_id = "xeos3111"
        slug = "random-slug"
        reason = "Model not found."

        # Setup mocks
        mock_instance = MagicMock()
        mock_instance.model_id = model_id
        mock_instance.slug = slug
        mock_model_base.return_value = mock_instance

        mock_fetch.return_value = FetchResult(False, reason)

        result = fetch(
            model=model_id,
            overwrite=True,
            from_dir=None,
            from_github=False,
            from_dockerhub=True,
            version="latest",
            from_s3=False,
            from_hosted=False,
            hosted_url=None,
            verbose=True,
        )

        logger.warning(f"❌ Fetch failed: {mock_fetch.return_value.reason}")

        # Assert the result signals failure
        assert result is False
        assert mock_fetch.return_value.fetch_success is False
        assert mock_fetch.return_value.reason == reason

        # Confirm the mock calls occurred
        mock_fetch.assert_called_once()
        mock_model_base.assert_called_once()


