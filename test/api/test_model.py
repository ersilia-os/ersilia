from unittest.mock import patch
from ersilia.api import Model

MODEL_ID = "eos3b5e"

INPUT_LIST = [
    "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O"
]

@patch('ersilia.api.create_api.Model._is_docker_running')
def test_fetch(mock_docker):
    """Test that fetch() succeeds when Docker is running"""
    mock_docker.return_value = True
    mdl = Model(MODEL_ID)
    result = mdl.fetch()
    assert result is True
    assert mdl.is_fetched()

@patch('ersilia.api.create_api.Model._is_docker_running')
def test_fetch_fails_when_docker_is_not_active(mock_docker):
    """Test that fetch() returns False when Docker is not running and from_dockerhub=True"""
    mock_docker.return_value = False
    mdl = Model(MODEL_ID)
    result = mdl.fetch()
    assert result is False

@patch('ersilia.api.create_api.Model._is_docker_running')
def test_serve_run_and_close(mock_docker):
    """Test that serve() returns correct dict and close() returns bool"""
    mock_docker.return_value = True
    mdl = Model(MODEL_ID)

    serve_result = mdl.serve()
    assert isinstance(serve_result, dict)
    assert "url" in serve_result
    assert "session" in serve_result
    assert "server" in serve_result

    df = mdl.run(input_list=INPUT_LIST)
    assert df.shape[0] == len(INPUT_LIST)

    close_result = mdl.close()
    assert isinstance(close_result, bool)

@patch('ersilia.api.commands.delete.delete')
def test_delete(mock_delete_cmd):
    """Test that delete() returns the result from underlying delete command"""
    mock_delete_cmd.return_value = True
    mdl = Model(MODEL_ID)
    result = mdl.delete()
    assert result is True
