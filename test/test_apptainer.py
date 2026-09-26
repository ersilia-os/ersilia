"""Fetching and serving models as Apptainer images, with Apptainer mocked."""

import asyncio
import os
from unittest.mock import MagicMock, patch

import pytest
import requests

from ersilia.setup.requirements.apptainer import ApptainerRequirement
from ersilia.utils import apptainer
from ersilia.utils.exceptions_utils.cli_exceptions import (
    ApptainerDownloadError,
    ApptainerImageNotFoundError,
    ApptainerNotInstalledError,
    ApptainerNotLinuxError,
    ModelStartError,
)

MODEL_ID = "eos4e40"


# Requirements


@patch("ersilia.setup.requirements.apptainer.platform.system", return_value="Darwin")
def test_not_linux_says_so(_system):
    with pytest.raises(ApptainerNotLinuxError) as e:
        ApptainerRequirement().check()
    assert "only be run on Linux, not on macOS" in e.value.message
    assert "DockerHub" in e.value.hints


@patch("ersilia.setup.requirements.apptainer.shutil.which", return_value=None)
@patch("ersilia.setup.requirements.apptainer.platform.system", return_value="Linux")
def test_not_installed_says_so(_system, _which):
    with pytest.raises(ApptainerNotInstalledError) as e:
        ApptainerRequirement().check()
    assert "apptainer.org" in e.value.hints


@patch("ersilia.setup.requirements.apptainer.platform.system", return_value="Linux")
def test_singularity_is_accepted(_system):
    which = {"singularity": "/usr/bin/singularity"}.get
    with patch("ersilia.setup.requirements.apptainer.shutil.which", side_effect=which):
        assert ApptainerRequirement().check() == "singularity"


# The image store


def test_urls_come_from_the_bucket_constant():
    with patch.object(apptainer, "SIF_BUCKET_URL", "https://example.org/sifs"):
        assert apptainer.sif_url(MODEL_ID, "v2") == (
            "https://example.org/sifs/eos4e40_v2.sif"
        )


def _head(status, length=None):
    response = MagicMock(status_code=status)
    response.headers = {"Content-Length": str(length)} if length else {}
    return response


def test_remote_size():
    with patch("ersilia.utils.apptainer.requests.head", return_value=_head(200, 42)):
        assert apptainer.remote_size(MODEL_ID, "v1") == 42
    # The bucket cannot be listed, so a missing image answers 403.
    with patch("ersilia.utils.apptainer.requests.head", return_value=_head(403)):
        assert apptainer.remote_size(MODEL_ID, "v9") is None


def test_unreachable_store_is_a_clear_error():
    error = requests.exceptions.ConnectionError("offline")
    with patch("ersilia.utils.apptainer.requests.head", side_effect=error):
        with pytest.raises(ApptainerDownloadError):
            apptainer.remote_size(MODEL_ID, "v1")


def test_latest_version_is_the_last_one_that_exists():
    sizes = {"v1": 10, "v2": 20}
    with patch.object(
        apptainer, "remote_size", side_effect=lambda m, v: sizes.get(v)
    ) as size:
        assert apptainer.latest_version(MODEL_ID) == "v2"
    # v1, v2, then v3 is missing: nothing further is asked.
    assert [c.args[1] for c in size.call_args_list] == ["v1", "v2", "v3"]


def test_model_without_image():
    with patch.object(apptainer, "remote_size", return_value=None):
        with pytest.raises(ApptainerImageNotFoundError) as e:
            apptainer.latest_version(MODEL_ID)
    assert "has no Apptainer image" in e.value.message


# Downloading


class _Download:
    """A streamed response that sends ``data`` in two chunks."""

    def __init__(self, data):
        self.data = data

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def raise_for_status(self):
        pass

    def iter_content(self, chunk_size):
        half = len(self.data) // 2
        yield self.data[:half]
        yield self.data[half:]


@pytest.fixture
def sif_dir(tmp_path):
    with patch.object(apptainer, "SIF_DIR", str(tmp_path)):
        yield tmp_path


def test_download_writes_the_image(sif_dir):
    data = b"x" * 1000
    with (
        patch.object(apptainer, "remote_size", return_value=len(data)),
        patch("ersilia.utils.apptainer.requests.get", return_value=_Download(data)),
    ):
        path = apptainer.download(MODEL_ID, "v1", verbose=True)
    assert path == str(sif_dir / "eos4e40_v1.sif")
    assert open(path, "rb").read() == data
    assert not os.path.exists(path + ".part")


def test_complete_image_is_not_downloaded_again(sif_dir):
    (sif_dir / "eos4e40_v1.sif").write_bytes(b"x" * 10)
    with (
        patch.object(apptainer, "remote_size", return_value=10),
        patch("ersilia.utils.apptainer.requests.get") as get,
    ):
        apptainer.download(MODEL_ID, "v1", verbose=True)
    get.assert_not_called()


def test_incomplete_download_leaves_nothing_behind(sif_dir):
    with (
        patch.object(apptainer, "remote_size", return_value=1000),
        patch("ersilia.utils.apptainer.requests.get", return_value=_Download(b"x" * 10)),
    ):
        with pytest.raises(ApptainerDownloadError):
            apptainer.download(MODEL_ID, "v1", verbose=True)
    assert list(sif_dir.iterdir()) == []


def test_download_of_a_missing_version():
    with patch.object(apptainer, "remote_size", return_value=None):
        with pytest.raises(ApptainerImageNotFoundError) as e:
            apptainer.download(MODEL_ID, "v7")
    assert "version v7" in e.value.message


# Running commands inside an image


def test_exec_args():
    assert apptainer.SimpleApptainer("apptainer").exec_args("/i.sif") == [
        "apptainer",
        "exec",
        "/i.sif",
    ]
    assert apptainer.SimpleApptainer("apptainer", use_unshare=True).exec_args(
        "/i.sif"
    ) == ["unshare", "-r", "apptainer", "exec", "/i.sif"]


# Fetching


def test_fetch_checks_requirements_before_downloading():
    from ersilia.hub.fetch.lazy_fetchers.apptainer import ModelApptainerFetcher

    with (
        patch(
            "ersilia.hub.fetch.lazy_fetchers.apptainer.ApptainerRequirement.check",
            side_effect=ApptainerNotLinuxError("macOS"),
        ),
        patch.object(apptainer, "download") as download,
        patch.object(apptainer, "latest_version") as latest,
    ):
        with pytest.raises(SystemExit):
            asyncio.run(ModelApptainerFetcher().fetch(MODEL_ID))
    download.assert_not_called()
    latest.assert_not_called()


# Serving


def _service(tmp_path, info):
    from ersilia.serve.services import ApptainerImageService

    service = ApptainerImageService.__new__(ApptainerImageService)
    service.model_id = MODEL_ID
    service.logger = MagicMock()
    service.info = info
    service.port = None
    service._port_given = False
    service.url = None
    service.pid = -1
    return service


def _info(tmp_path):
    sif = tmp_path / "eos4e40_v1.sif"
    sif.write_bytes(b"sif")
    return {
        "apptainer": True,
        "sif_path": str(sif),
        "bundle_path": "/opt/ersilia/bundles/eos4e40",
        "binary": "apptainer",
        "use_unshare": False,
    }


@patch("ersilia.setup.requirements.apptainer.ApptainerRequirement.check")
def test_serve_records_the_server_not_the_wrapper(_check, tmp_path):
    service = _service(tmp_path, _info(tmp_path))
    wrapper = MagicMock(pid=111)
    wrapper.poll.return_value = 0  # 'apptainer' exits once the server is up
    with (
        patch("subprocess.Popen", return_value=wrapper) as popen,
        patch("ersilia.serve.services.find_free_port", return_value=8123),
        patch("ersilia.serve.services.make_temp_dir", return_value=str(tmp_path)),
        patch.object(service, "_is_ready", return_value=True),
        patch.object(service, "_listening_pid", return_value=222),
        patch.object(service, "_get_apis_from_apis_list", return_value=None),
    ):
        service.serve()
    command = popen.call_args.args[0]
    assert command[:3] == ["apptainer", "exec", _info(tmp_path)["sif_path"]]
    assert command[3:] == [
        "ersilia_model_serve",
        "--bundle_path",
        "/opt/ersilia/bundles/eos4e40",
        "--port",
        "8123",
    ]
    assert service.pid == 222
    assert service.url == "http://127.0.0.1:8123"
    assert service._apis_list == ["run"]


@patch("ersilia.setup.requirements.apptainer.ApptainerRequirement.check")
def test_serve_reports_a_crash_and_cleans_up(_check, tmp_path):
    service = _service(tmp_path, _info(tmp_path))
    (tmp_path / "serve.log").write_text("")

    def crash(*args, **kwargs):
        (tmp_path / "serve.log").write_text(
            "Traceback (most recent call last):\n  ...\nImportError: boom\n"
        )
        return MagicMock(pid=111, poll=MagicMock(return_value=None))

    with (
        patch("subprocess.Popen", side_effect=crash),
        patch("ersilia.serve.services.find_free_port", return_value=8123),
        patch("ersilia.serve.services.make_temp_dir", return_value=str(tmp_path)),
        patch.object(service, "_is_ready", return_value=False),
        patch.object(service, "_stop") as stop,
    ):
        with pytest.raises(ModelStartError) as e:
            service.serve()
    assert "stopped while starting" in e.value.message
    stop.assert_called_once()


@patch("ersilia.setup.requirements.apptainer.ApptainerRequirement.check")
def test_serve_without_the_image(_check, tmp_path):
    info = dict(_info(tmp_path), sif_path=str(tmp_path / "gone.sif"))
    service = _service(tmp_path, info)
    with pytest.raises(ModelStartError) as e:
        service.serve()
    assert "--from_apptainer" in e.value.hints


def test_not_available_unless_fetched_as_apptainer(tmp_path):
    assert not _service(tmp_path, None).is_available()
    with patch(
        "ersilia.setup.requirements.apptainer.shutil.which", return_value="/x"
    ):
        assert _service(tmp_path, _info(tmp_path)).is_available()
