import pytest
from unittest.mock import MagicMock, patch

# ← fix import to point at the commands sub‑package
from ersilia.api.commands.serve import serve

MODEL = "eos3b5e"
URL = "http://localhost"
SESSION_DIR = "/tmp/session"
SCL = "pulled_docker"
SLUG = "molecular-weight"


@patch("ersilia.api.commands.serve.logger")
@patch("ersilia.api.commands.serve.SetupRedis")
@patch("ersilia.api.commands.serve.register_model_session")
@patch("ersilia.api.commands.serve.echo")
@patch("ersilia.api.commands.serve.ErsiliaModel")
def test_serve_success(
    mock_model_cls, mock_echo, mock_register, mock_setupredis, mock_logger
):
    # Arrange: a valid model instance
    inst = MagicMock()
    inst.is_valid.return_value = True
    inst.url = URL
    inst.session._session_dir = SESSION_DIR
    inst.scl = SCL
    inst.model_id = MODEL
    inst.slug = SLUG
    inst.get_apis.return_value = ["run"]

    mock_model_cls.return_value = inst

    # Act
    result = serve(MODEL)

    # Assert
    mock_logger.set_verbosity.assert_called_once_with(0)
    mock_model_cls.assert_called_once_with(
        MODEL,
        output_source=None,
        preferred_port=None,
        cache=True,
        maxmemory=None,
    )
    mock_setupredis.assert_called_once_with(True, None)
    inst.is_valid.assert_called_once()
    inst.serve.assert_called_once_with(track_runs=None)
    mock_register.assert_called_once_with(MODEL, SESSION_DIR)
    assert result == (URL, SESSION_DIR, SCL)


def test_serve_invalid_memory_fraction():
    # max_cache_memory_frac outside [0.2,0.7] should error
    with pytest.raises(RuntimeError) as exc:
        serve(MODEL, max_cache_memory_frac=0.1)
    assert "outside of recommended range" in str(exc.value)


@patch("ersilia.api.commands.serve.ErsiliaModel")
def test_serve_invalid_model(mock_model_cls):
    # Arrange: model.is_valid() -> False
    inst = MagicMock()
    inst.is_valid.return_value = False
    mock_model_cls.return_value = inst

    with pytest.raises(RuntimeError) as exc:
        serve(MODEL)
    assert f"Model {inst.model_id} is not valid" in str(exc.value)


@patch("ersilia.api.commands.serve.ErsiliaModel")
def test_serve_no_url(mock_model_cls):
    # Arrange: valid but no URL
    inst = MagicMock()
    inst.is_valid.return_value = True
    inst.url = None
    mock_model_cls.return_value = inst

    with pytest.raises(RuntimeError) as exc:
        serve(MODEL)
    assert "No URL found. Service unsuccessful." in str(exc.value)
