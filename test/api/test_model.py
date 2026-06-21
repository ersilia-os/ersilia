from unittest.mock import patch
import pandas as pd
from ersilia.api import Model

MODEL_ID = "eos3b5e"

INPUT_LIST = [
    "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O"
]

@patch('ersilia.api.create_api.Model._is_docker_running')
@patch('ersilia.api.commands.fetch.fetch')
def test_fetch(mock_fetch_cmd, mock_docker):
    """Test that fetch() returns result from underlying fetch command"""
    mock_docker.return_value = True
    mock_fetch_cmd.return_value = True
    mdl = Model(MODEL_ID)
    result = mdl.fetch()
    assert result is True

@patch('ersilia.api.create_api.Model._is_docker_running')
def test_fetch_fails_when_docker_is_not_active(mock_docker):
    """Test that fetch() returns False when Docker is not running and from_dockerhub=True"""
    mock_docker.return_value = False
    mdl = Model(MODEL_ID)
    result = mdl.fetch()
    assert result is False

@patch('ersilia.api.commands.serve.serve')
def test_serve_returns_dict(mock_serve_cmd):
    """Test that serve() returns correct dict with url, session, and server keys"""
    mock_serve_cmd.return_value = ("http://localhost:5000", "test_session", "test_server")
    mdl = Model(MODEL_ID)
    serve_result = mdl.serve()

    assert isinstance(serve_result, dict)
    assert "url" in serve_result
    assert "session" in serve_result
    assert "server" in serve_result
    assert serve_result["url"] == "http://localhost:5000"

@patch('ersilia.api.commands.run.run')
def test_run_returns_dataframe(mock_run_cmd):
    """Test that run() returns a pandas dataframe"""
    mock_df = pd.DataFrame({"prediction": [0.5]})
    mock_run_cmd.return_value = mock_df

    mdl = Model(MODEL_ID)
    result = mdl.run(input_list=INPUT_LIST)

    assert isinstance(result, pd.DataFrame)
    assert result.shape[0] == 1

@patch('ersilia.api.commands.close.close')
def test_close(mock_close_cmd):
    """Test that close() returns the result from underlying close command"""
    mock_close_cmd.return_value = True
    mdl = Model(MODEL_ID)
    close_result = mdl.close()

    assert isinstance(close_result, bool)
    assert close_result is True

@patch('ersilia.api.commands.delete.delete')
def test_delete(mock_delete_cmd):
    """Test that delete() returns the result from underlying delete command"""
    mock_delete_cmd.return_value = True
    mdl = Model(MODEL_ID)
    result = mdl.delete()
    assert result is True
