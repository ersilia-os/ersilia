import pytest
import random
from unittest.mock import patch, Mock, AsyncMock, PropertyMock
from click.testing import CliRunner
from ersilia.cli.commands.run import run_cmd
from ersilia.serve.standard_api import StandardCSVRunApi
from ersilia.serve.autoservice import AutoService
from ersilia.serve.services import DockerImageService
from ersilia.core.session import Session
from ersilia.core.model import ErsiliaModel
from ersilia.utils.logging import logger
from ersilia.hub.fetch.fetch import ModelFetcher
from .utils import create_compound_input_csv

URL = "http://localhost"
PORT = 8001
MODEL_ID = "eos3b5e"
API_NAME = "run"
INPUT = "NCCCCCCCCCCNS(=O)(=O)c1cccc2c(Cl)cccc12"
INPUT_CSV = "input.csv"
RESULT_CSV = "result.csv"
MIN_WEIGHT = 40.0
MAX_WEIGHT = 60.0
HEADER = ["key", "input", "value"]


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
def mock_convn_api_get_apis():
    def mock_get_api_side_effect():
        return [API_NAME]

    with patch.object(
        ErsiliaModel, "get_apis", side_effect=mock_get_api_side_effect
    ) as mock_get_apis:
        yield mock_get_apis

@pytest.fixture
def mock_get_url():
    with patch.object(
        ErsiliaModel, 
        "_get_url", 
        return_value=f"{URL}:{PORT}"
        ) as mock_url:
        yield mock_url

@pytest.fixture
def mock_get_input():
    with patch.object(
        StandardCSVRunApi, 
        "get_input_type",
          return_value=["Compound"]
    ) as mock_input_type:
        yield mock_input_type


@pytest.fixture
def mock_std_header():
    with patch.object(
        StandardCSVRunApi, 
        "get_expected_output_header", 
        return_value=HEADER
    ) as mock_method:
        yield mock_method


@pytest.fixture
def mock_is_ready():
    with patch.object(
        StandardCSVRunApi, 
        "is_ready", 
        return_value=True
    ) as mock_is_ready:
        yield mock_is_ready


@pytest.fixture
def mock_is_amenable():
    with patch.object(
        StandardCSVRunApi, 
        "is_amenable", 
        return_value=True
    ) as mock_is_amenable:
        yield mock_is_amenable

@pytest.fixture
def compound_csv():
    create_compound_input_csv(INPUT_CSV)
    yield INPUT_CSV

@pytest.fixture
def mock_std_api_post(compound_csv):
    def mock_post_side_effect(input, output, output_source):
        api_instance = StandardCSVRunApi(
            model_id=MODEL_ID, 
            url=URL
        )
        logger.info(f"Input: {input}")
        input_data = api_instance.serialize_to_json(input)

        logger.info(f"Serialized Input Data: {input_data}")

        return_value = [
            {"value": round(random.uniform(MIN_WEIGHT, MAX_WEIGHT), 4)} 
            for _ in input_data
        ]

        try:
            csv_result = api_instance.serialize_to_csv(
                input_data, 
                return_value, output
            )
            logger.info(f"CSV Result: {csv_result}")
            return csv_result
        except Exception as e:
            logger.error(f"Exception in serialize_to_csv: {e}")
            return "Error"

    with patch.object(
        StandardCSVRunApi, "post", side_effect=mock_post_side_effect
    ) as mock_post:
        yield mock_post


@pytest.fixture
def mock_session():
    with patch.object(
        Session, 
        "current_model_id", 
        return_value=MODEL_ID
        ), \
    patch.object(
        Session,
        "current_service_class", 
        return_value="pulled_docker" 
    ), \
    patch.object(
        Session, 
        "tracking_status", 
        return_value=False), \
    patch.object(
        Session, 
        "current_output_source", 
        return_value="LOCAL_ONLY"
    ):
        yield

# For Standard API
def test_standard_api_string(
    mock_std_api_post,
    mock_fetcher,
    mock_set_apis,
    mock_get_input,
    mock_is_ready,
    mock_std_header,
    mock_is_amenable,
    mock_get_url,
    mock_session,
):
    runner = CliRunner()

    input_arg = INPUT      
    output_arg = RESULT_CSV
    result = runner.invoke(
        run_cmd(), ["-i", input_arg, "-o", output_arg]
    )
    assert result.exit_code == 0
    assert mock_get_url.called
    assert mock_get_input.called
    assert mock_std_header.called
    assert mock_is_ready.called
    assert mock_is_amenable.called
    assert mock_std_api_post.called
    assert mock_set_apis.called
    assert RESULT_CSV in result.output


def test_standard_api_csv(
    mock_std_api_post,
    mock_fetcher,
    mock_set_apis,
    mock_get_input,
    mock_is_ready,
    mock_std_header,
    mock_is_amenable,
    mock_get_url,
    mock_session
):
    runner = CliRunner()
    input_arg = INPUT_CSV
    output_arg = RESULT_CSV
    result = runner.invoke(
        run_cmd(), ["-i", input_arg, "-o", output_arg]
    )

    assert result.exit_code == 0
    assert mock_get_input.called
    assert mock_get_url.called
    assert mock_std_header.called
    assert mock_is_ready.called
    assert mock_is_amenable.called
    assert mock_std_api_post.called
    assert mock_set_apis.called
    assert RESULT_CSV in result.output

# For Conventional Run
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
def mock_set_apis():
    with patch.object(ErsiliaModel, "_set_apis", return_value=None) as mock_set_apis:
        yield mock_set_apis


def test_conv_api_string(
    mock_convn_api_get_apis, mock_api_task, mock_set_apis, mock_fetcher, mock_session
):
    runner = CliRunner()

    input_arg = INPUT
    output_arg = RESULT_CSV
    batch_size = 10
    result = runner.invoke(
        run_cmd(), ["-i", input_arg, "-b", str(batch_size)]
    )

    assert result.exit_code == 0
    assert mock_convn_api_get_apis.called
    assert mock_set_apis.called

def test_conv_api_csv(
    mock_convn_api_get_apis, 
    mock_api_task,
    mock_set_apis,
    mock_fetcher,
    mock_session
):
    runner = CliRunner()
    input_arg = INPUT_CSV
    output_arg = RESULT_CSV
    batch_size = 10
    result = runner.invoke(
        run_cmd(), ["-i", input_arg, "-b", str(batch_size)]
    )
    logger.info(result.output)
    assert result.exit_code == 0
    assert mock_convn_api_get_apis.called
    assert mock_set_apis.called
