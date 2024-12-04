import asyncio
import pytest
import random
from unittest.mock import patch, AsyncMock
from ersilia.hub.fetch.fetch import ModelFetcher
from ersilia import ErsiliaModel
from ersilia.core.session import Session

MODELS = ["eos0t01", "eos3b5e", "eos0t03", "eos0t04"]
RESULTS = [0, 312.89, 0, 0]
API_NAME="run"
MIN_WEIGHT = 40.0
MAX_WEIGHT = 60.0

@pytest.fixture
def mock_fetcher():
    with patch(
        "ersilia.hub.fetch.fetch.ModelFetcher.fetch", 
        new_callable=AsyncMock
    ) as mock_fetch:
        yield mock_fetch

@pytest.fixture
def mock_set_apis():
    with patch.object(
        ErsiliaModel, 
        "_set_apis", 
        return_value=None
    ) as mock_set_apis:
        yield mock_set_apis

@pytest.fixture
def mock_session():
    with patch.object(Session, "current_model_id", return_value=MODELS[1]), \
         patch.object(Session, "current_service_class", return_value="docker"), \
         patch.object(Session, "tracking_status", return_value=False), \
         patch.object(Session, "current_output_source", return_value="LOCAL_ONLY"):
        yield

@pytest.fixture
def mock_convn_api_get_apis():
    def mock_get_api_side_effect():
        return [API_NAME]

    with patch.object(
        ErsiliaModel, "get_apis", side_effect=mock_get_api_side_effect
    ) as mock_get_apis:
        yield mock_get_apis


@pytest.fixture
def mock_api_task():
    def mock_api_task_side_effect(api_name, input, output, batch_size):
        if output is None:
            length = len(input)

            def result_generator():
                for _ in range(length):
                    yield {"value": round(random.uniform(MIN_WEIGHT, MAX_WEIGHT), 3)}

            return result_generator()
        return output

    with patch.object(
        ErsiliaModel, "api_task", side_effect=mock_api_task_side_effect
    ) as mock_task:
        yield mock_task

@pytest.fixture
def mock_close():
    with patch.object(
        ErsiliaModel, 
        "close", 
        return_value=None
        ) as mock_close_:
        yield mock_close_

@pytest.fixture
def mock_serve():
    with patch.object(
        ErsiliaModel, 
        "serve", 
        return_value=None
        ) as mock_serve_:
        yield mock_serve_

@pytest.fixture
def mock_run():
    with patch.object(
        ErsiliaModel, 
        "run", 
        return_value=RESULTS[1]
        ) as mock_run_:
        yield mock_run_

@patch("ersilia.core.model.ErsiliaModel")
def test_models(
    mock_ersilia_model, 
    mock_fetcher, 
    mock_session,
    mock_set_apis,
    mock_convn_api_get_apis,
    mock_api_task,
    mock_serve,
    mock_run,
    mock_close
):
    MODEL_ID = MODELS[1]
    INPUT = "CCCC"

    mf = ModelFetcher(overwrite=True)
    asyncio.run(mf.fetch(MODEL_ID))
    
    em = ErsiliaModel(
        model=MODEL_ID,
        service_class="docker",
        output_source="LOCAL_ONLY")
    
    result = em.run(
        input=INPUT,
        output="result.csv",
        batch_size=100,
        track_run=False,
        try_standard=False,
    )

    em.serve()
    em.close()

    assert result == RESULTS[1]
    assert mock_fetcher.called
    assert mock_serve.called
    assert mock_run.called
    assert mock_close.called
