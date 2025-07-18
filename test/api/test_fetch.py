from unittest.mock import MagicMock, patch
import pytest

from ersilia.api.commands.fetch import fetch  # API fetch function
from ersilia.hub.fetch.fetch import FetchResult
from ersilia.utils.logging import logger

MODEL_ID = "eos3b5e"
SLUG = "molecular-weight"


# ✅ Fixture to mock ModelFetcher.fetch (async method)
@pytest.fixture
def mock_model_fetcher():
    with patch("ersilia.api.commands.fetch.ModelFetcher.fetch") as mock_fetch:
        mock_fetch.return_value = FetchResult(True, "Model fetched successfully.")
        yield mock_fetch


# ✅ Fixture to mock ModelBase
@pytest.fixture
def mock_model_base():
    with patch("ersilia.api.commands.fetch.ModelBase") as MockModelBase:
        mock_instance = MagicMock()
        mock_instance.model_id = MODEL_ID
        mock_instance.slug = SLUG
        MockModelBase.return_value = mock_instance
        yield MockModelBase


# ✅ Successful fetch test using API
def test_fetch_api_success(mock_model_fetcher, mock_model_base):
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

    logger.info(f"✅ Fetch succeeded: {mock_model_fetcher.return_value.reason}")

    assert result is not False
    mock_model_fetcher.assert_called_once()
    mock_model_base.assert_called_once()


# ❌ Failure fetch test with SystemExit handling
def test_fetch_api_failure():
    with patch("ersilia.api.commands.fetch.ModelBase") as mock_model_base, \
         patch("ersilia.api.commands.fetch.ModelFetcher.fetch") as mock_fetch:

        model_id = "xeos3111"
        slug = "random-slug"
        reason = "Model not found."

        # Simulate invalid model
        mock_instance = MagicMock()
        mock_instance.model_id = model_id
        mock_instance.slug = slug
        mock_model_base.return_value = mock_instance
        mock_fetch.return_value = FetchResult(False, reason)

        
        fetch(
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
        mock_fetch.assert_called_once()
        mock_model_base.assert_called_once()


