import os
from unittest.mock import MagicMock, call, patch

from ersilia.utils.docker import ERSILIA_USER_LABEL, SimpleDocker


def _mock_subprocess_result(stdout="", stderr=b"", returncode=0):
    result = MagicMock()
    result.stdout = stdout
    result.stderr = stderr
    result.returncode = returncode
    return result


CURRENT_UID = str(os.getuid())
OTHER_UID = str(os.getuid() + 1)


class TestIsOwnedByCurrentUser:
    def test_label_matches_current_user(self):
        with patch("subprocess.run", return_value=_mock_subprocess_result(stdout=CURRENT_UID)):
            d = SimpleDocker()
            assert d.is_owned_by_current_user("ersiliaos", "eos0001", "latest") is True

    def test_label_belongs_to_different_user(self):
        with patch("subprocess.run", return_value=_mock_subprocess_result(stdout=OTHER_UID)):
            d = SimpleDocker()
            assert d.is_owned_by_current_user("ersiliaos", "eos0001", "latest") is False

    def test_no_label_returns_true_for_backward_compat(self):
        with patch("subprocess.run", return_value=_mock_subprocess_result(stdout="")):
            d = SimpleDocker()
            assert d.is_owned_by_current_user("ersiliaos", "eos0001", "latest") is True


class TestLabelWithCurrentUser:
    def test_calls_docker_build_with_correct_label(self):
        inspect_result = _mock_subprocess_result(stdout="sha256:abc123")
        build_result = _mock_subprocess_result()
        rmi_result = _mock_subprocess_result()

        with patch("subprocess.run", side_effect=[inspect_result, build_result, rmi_result]) as mock_run:
            d = SimpleDocker()
            d.label_with_current_user("ersiliaos", "eos0001", "latest")

            calls = mock_run.call_args_list
            build_args = calls[1][0][0]
            assert "build" in build_args
            assert f"{ERSILIA_USER_LABEL}={CURRENT_UID}" in " ".join(build_args)
            assert "ersiliaos/eos0001:latest" in build_args

    def test_removes_dangling_image_after_labeling(self):
        old_id = "sha256:abc123"
        inspect_result = _mock_subprocess_result(stdout=old_id)
        build_result = _mock_subprocess_result()
        rmi_result = _mock_subprocess_result()

        with patch("subprocess.run", side_effect=[inspect_result, build_result, rmi_result]) as mock_run:
            d = SimpleDocker()
            d.label_with_current_user("ersiliaos", "eos0001", "latest")

            rmi_call_args = mock_run.call_args_list[2][0][0]
            assert "rmi" in rmi_call_args
            assert old_id in rmi_call_args

    def test_logs_warning_on_failure(self):
        inspect_result = _mock_subprocess_result(stdout="sha256:abc123")
        build_result = _mock_subprocess_result(stderr=b"error message", returncode=1)

        with patch("subprocess.run", side_effect=[inspect_result, build_result]):
            d = SimpleDocker()
            with patch.object(d.logger, "warning") as mock_warn:
                d.label_with_current_user("ersiliaos", "eos0001", "latest")
                mock_warn.assert_called_once()
                assert "Could not label" in mock_warn.call_args[0][0]


class TestDelete:
    def test_skips_deletion_for_other_users_image(self):
        with patch("subprocess.run", return_value=_mock_subprocess_result(stdout=OTHER_UID)):
            d = SimpleDocker()
            with patch("ersilia.utils.docker.run_command") as mock_run_cmd:
                with patch.object(d.logger, "warning") as mock_warn:
                    d.delete("ersiliaos", "eos0001", "latest")
                    mock_run_cmd.assert_not_called()
                    mock_warn.assert_called_once()
                    assert "different user" in mock_warn.call_args[0][0]

    def test_deletes_own_image(self):
        with patch("subprocess.run", return_value=_mock_subprocess_result(stdout=CURRENT_UID)):
            d = SimpleDocker()
            with patch("ersilia.utils.docker.run_command") as mock_run_cmd:
                d.delete("ersiliaos", "eos0001", "latest")
                mock_run_cmd.assert_called_once()
                assert "docker rmi" in mock_run_cmd.call_args[0][0]

    def test_deletes_unlabeled_image(self):
        with patch("subprocess.run", return_value=_mock_subprocess_result(stdout="")):
            d = SimpleDocker()
            with patch("ersilia.utils.docker.run_command") as mock_run_cmd:
                d.delete("ersiliaos", "eos0001", "latest")
                mock_run_cmd.assert_called_once()


class TestPulledDockerImageServiceEnv:
    """The container does not inherit the shell's environment, so serve() must
    forward the variables the model server reads."""

    def _serve_and_get_env(self, monkeypatch):
        from ersilia.serve.services import PulledDockerImageService

        svc = PulledDockerImageService.__new__(PulledDockerImageService)
        svc.model_id = "eos0000"
        svc.image_name = "ersiliaos/eos0000:latest"
        svc.port = 8080
        svc._port_given = False
        svc._mem_gb = None
        svc.logger = MagicMock()
        svc.client = MagicMock()
        monkeypatch.setattr(svc, "_create_docker_network", lambda: None)
        monkeypatch.setattr(svc, "_wait_until_container_is_running", lambda: None)
        monkeypatch.setattr(svc, "_get_apis", lambda: [])
        svc.serve()
        return svc.client.containers.run.call_args.kwargs["environment"]

    def test_forwards_ersilia_warm_when_set(self, monkeypatch):
        monkeypatch.setenv("ERSILIA_WARM", "0")
        env = self._serve_and_get_env(monkeypatch)
        assert env["ERSILIA_WARM"] == "0"

    def test_omits_ersilia_warm_when_unset(self, monkeypatch):
        monkeypatch.delenv("ERSILIA_WARM", raising=False)
        env = self._serve_and_get_env(monkeypatch)
        assert "ERSILIA_WARM" not in env
