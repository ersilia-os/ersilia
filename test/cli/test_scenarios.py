"""Scenarios where the CLI must guide the user; no Docker or network needed."""

from unittest.mock import MagicMock, patch

import pytest

from ersilia.serve.services import PulledDockerImageService
from ersilia.utils.exceptions_utils.cli_exceptions import (
    ModelStartError,
    PortInUseError,
)


def _service(status, state=None):
    svc = PulledDockerImageService.__new__(PulledDockerImageService)
    svc.model_id = "eos3b5e"
    svc.url = "http://0.0.0.0:1"
    svc.logger = MagicMock()
    svc.container = MagicMock(status=status, attrs={"State": state or {}})
    svc.container.logs.return_value = b"model log line"
    svc.is_url_available = lambda url: False
    return svc


def test_start_reports_a_container_killed_for_memory():
    svc = _service("exited", {"OOMKilled": True, "ExitCode": 137})
    with pytest.raises(ModelStartError, match="out of memory") as e:
        svc._wait_until_container_is_running()
    assert "Give Docker more memory" in e.value.hints


def test_start_reports_a_crashed_container():
    svc = _service("exited", {"OOMKilled": False, "ExitCode": 132})
    with pytest.raises(ModelStartError, match="exit code 132") as e:
        svc._wait_until_container_is_running()
    assert "ersilia -v serve eos3b5e" in e.value.hints


def test_start_gives_up_after_the_deadline():
    svc = _service("running")
    with patch("time.sleep"), pytest.raises(ModelStartError, match="did not start"):
        svc._wait_until_container_is_running(timeout=0)


def test_start_returns_once_the_server_answers():
    svc = _service("running")
    svc.is_url_available = lambda url: True
    svc._wait_until_container_is_running()


def test_a_port_given_by_the_user_must_be_free():
    svc = PulledDockerImageService.__new__(PulledDockerImageService)
    svc.model_id, svc.port, svc._port_given = "eos3b5e", 8080, True
    svc.image_name, svc._mem_gb = "ersiliaos/eos3b5e:latest", None
    svc.logger, svc.client = MagicMock(), MagicMock()
    svc._create_docker_network = lambda: None
    with (
        patch("ersilia.utils.ports.is_port_in_use", return_value=True),
        pytest.raises(PortInUseError, match="Port 8080 is already in use"),
    ):
        svc.serve()
    svc.client.containers.run.assert_not_called()


def test_a_failed_start_removes_its_container():
    svc = PulledDockerImageService.__new__(PulledDockerImageService)
    svc.model_id, svc.port, svc._port_given = "eos3b5e", 8080, False
    svc.image_name, svc._mem_gb = "ersiliaos/eos3b5e:latest", None
    svc.logger, svc.client = MagicMock(), MagicMock()
    svc._create_docker_network = lambda: None
    svc._wait_until_container_is_running = MagicMock(side_effect=KeyboardInterrupt)
    with (
        patch("ersilia.serve.services.stop_containers_by_name") as stop,
        pytest.raises(KeyboardInterrupt),
    ):
        svc.serve()
    stop.assert_called_once_with([svc.container_name])
