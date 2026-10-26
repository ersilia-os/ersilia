from unittest.mock import MagicMock, patch

import docker
import pytest
from click.testing import CliRunner

from ersilia.cli.commands.fetch import fetch_cmd
from ersilia.hub.content.information import version_label
from ersilia.hub.fetch.fetch import FetchResult
from ersilia.publish.test.services.constants import STATUS_CONFIGS
from ersilia.publish.test.services.inspect import ModelInspector
from ersilia.serve.services import PulledDockerImageService
from ersilia.utils.docker import gpu_image_tag, is_gpu_image_tag
from ersilia.utils.exceptions_utils.cli_exceptions import ModelStartError

MODEL_ID = "eos9o72"


@pytest.mark.parametrize(
    "tag, expected",
    [
        (None, "latest-gpu"),
        ("latest", "latest-gpu"),
        ("v1.0.0", "v1.0.0-gpu"),
        ("latest-gpu", "latest-gpu"),
    ],
)
def test_gpu_image_tag(tag, expected):
    assert gpu_image_tag(tag) == expected
    assert is_gpu_image_tag(expected)


def test_cpu_tags_are_not_gpu():
    assert not is_gpu_image_tag("latest")
    assert not is_gpu_image_tag(None)


def test_version_label_gpu():
    card = {"Release": "v1.0.1"}
    assert version_label("latest-gpu", card) == "latest-gpu (v1.0.1)"
    assert version_label("latest", card) == "latest (v1.0.1)"
    assert version_label("v1.0.0-gpu", card) == "v1.0.0-gpu"


def _fetch(flags):
    mdl = MagicMock(model_id=MODEL_ID, slug="chemeleon")
    with (
        patch("ersilia.ModelBase", return_value=mdl),
        patch("ersilia.hub.fetch.fetch.ModelFetcher") as fetcher,
    ):
        fetcher.return_value.fetch.return_value = FetchResult(
            True, "Model fetched successfully"
        )
        with patch("nest_asyncio.apply"), patch("asyncio.run", side_effect=lambda r: r):
            result = CliRunner().invoke(fetch_cmd(), [MODEL_ID] + flags)
    return result, fetcher


@pytest.mark.parametrize(
    "flags, img_version",
    [
        ([], None),
        (["--gpu"], "latest-gpu"),
        (["--gpu", "--version", "v1.0.0"], "v1.0.0-gpu"),
    ],
)
def test_fetch_gpu_tag(flags, img_version):
    result, fetcher = _fetch(flags)
    assert result.exit_code == 0, result.output
    assert fetcher.call_args.kwargs["img_version"] == img_version


def test_fetch_gpu_needs_dockerhub():
    result, fetcher = _fetch(["--gpu", "--from_github"])
    assert result.exit_code == 2
    assert "only applies to DockerHub" in " ".join(result.output.split())
    fetcher.assert_not_called()


def _service(tag):
    srv = PulledDockerImageService.__new__(PulledDockerImageService)
    srv.model_id = MODEL_ID
    srv.docker_tag = tag
    srv.image_name = "ersiliaos/{0}:{1}".format(MODEL_ID, tag)
    srv.port = 8000
    srv._port_given = False
    srv._mem_gb = None
    srv.client = MagicMock()
    srv.logger = MagicMock()
    srv._create_docker_network = MagicMock()
    srv._wait_until_container_is_running = MagicMock()
    srv._get_apis = MagicMock(return_value=["run"])
    return srv


def test_serve_gpu_image_requests_gpus():
    srv = _service("latest-gpu")
    srv.serve()
    requests = srv.client.containers.run.call_args.kwargs["device_requests"]
    assert requests[0]["Count"] == -1
    assert requests[0]["Capabilities"] == [["gpu"]]


def test_serve_cpu_image_does_not_request_gpus():
    srv = _service("latest")
    srv.serve()
    assert "device_requests" not in srv.client.containers.run.call_args.kwargs


@pytest.mark.parametrize(
    "error",
    [
        'could not select device driver "" with capabilities: [[gpu]]',
        "failed to discover GPU vendor from CDI: no known GPU vendor found",
    ],
)
def test_serve_gpu_image_without_gpu_support(error):
    srv = _service("latest-gpu")
    srv.client.containers.run.side_effect = docker.errors.APIError(error)
    with patch("ersilia.serve.services.stop_containers_by_name"):
        with pytest.raises(ModelStartError, match="needs a GPU"):
            srv.serve()


CPU_COMMANDS = """commands:
    - ["pip", "torch", "2.10.0", "--index-url", "https://download.pytorch.org/whl/cpu"]
    - ["pip", "chemprop", "2.2.0"]
"""
GPU_COMMANDS = """gpu:
    - ["pip", "torch", "2.10.0", "--index-url", "https://download.pytorch.org/whl/cu126"]
    - ["pip", "chemprop", "2.2.0"]
"""
GPU_RUN_SH = (
    "MAIN=main.py\n"
    'if [ "$ERSILIA_DEVICE" = "gpu" ]; then MAIN=main_gpu.py; fi\n'
    "python $1/code/$MAIN $2 $3\n"
)


def _model_dir(tmp_path, gpu=GPU_COMMANDS, run_sh=GPU_RUN_SH, main_gpu=True, device=None):
    code = tmp_path / "model" / "framework" / "code"
    code.mkdir(parents=True)
    (code / "main.py").write_text("")
    if main_gpu:
        (code / "main_gpu.py").write_text("")
    (tmp_path / "model" / "framework" / "run.sh").write_text(run_sh)
    (tmp_path / "install.yml").write_text('python: "3.11"\n' + CPU_COMMANDS + (gpu or ""))
    (tmp_path / "metadata.yml").write_text(f"Device: {device}\n" if device else "")
    return ModelInspector(MODEL_ID, str(tmp_path))


def test_gpu_build_complete(tmp_path):
    inspector = _model_dir(tmp_path, device="[CPU, GPU]")
    assert inspector.check_gpu_build().success is True
    assert inspector.check_dependencies_are_valid().success


def test_gpu_build_skipped_without_gpu_block(tmp_path):
    result = _model_dir(tmp_path, gpu=None).check_gpu_build()
    assert result.success == str(STATUS_CONFIGS.SKIPPED)


def test_gpu_device_without_gpu_block(tmp_path):
    assert _model_dir(tmp_path, gpu=None, device="[GPU]").check_gpu_build().success is False


def test_gpu_build_needs_main_gpu(tmp_path):
    result = _model_dir(tmp_path, main_gpu=False).check_gpu_build()
    assert result.success is False and "main_gpu.py is missing" in result.details


def test_gpu_build_needs_run_sh_switch(tmp_path):
    result = _model_dir(tmp_path, run_sh="python $1/code/main.py $2 $3\n").check_gpu_build()
    assert result.success is False and "ERSILIA_DEVICE" in result.details


def test_gpu_build_warns_on_drift(tmp_path):
    result = _model_dir(tmp_path, gpu=GPU_COMMANDS.replace("2.2.0", "2.3.0")).check_gpu_build()
    assert result.success == str(STATUS_CONFIGS.WARNING)
    assert "chemprop" in result.details and "torch" not in result.details


@pytest.mark.parametrize(
    "gpu, error",
    [
        ("gpu: []\n", "non-empty list"),
        ('gpu:\n    - ["pip", "rdkit"]\n', "gpu: Missing version"),
    ],
)
def test_gpu_block_format(tmp_path, gpu, error):
    result = _model_dir(tmp_path, gpu=gpu).check_dependencies_are_valid()
    assert result.success is False and error in result.details
