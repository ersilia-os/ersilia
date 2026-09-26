"""Unit tests of the Python API. They mock the core, so no Docker or network."""

import glob
import inspect
import os
import subprocess
import sys
import tempfile
from collections import namedtuple
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from ersilia.api import Catalog, ErsiliaError, Model
from ersilia.api._runtime import library_call
from ersilia.cli.commands.catalog import catalog_cmd
from ersilia.cli.commands.example import example_cmd
from ersilia.cli.commands.fetch import fetch_cmd
from ersilia.cli.commands.info import info_cmd
from ersilia.cli.commands.run import run_cmd
from ersilia.cli.commands.serve import serve_cmd
from ersilia.utils import terminal
from ersilia.utils.echo import echo
from ersilia.utils.exceptions_utils.api_exceptions import (
    InvalidOptionError,
    ModelFetchError,
    ModelNotServedError,
    SessionBusyError,
)
from ersilia.utils.exceptions_utils.exceptions import InvalidModelIdentifierError
from ersilia.utils.exceptions_utils.throw_ersilia_exception import (
    throw_ersilia_exception,
)

FetchResult = namedtuple("FetchResult", ["fetch_success", "reason"])


def _model(model_id="eos3b5e"):
    # A Model without resolving the identifier (which may use the network).
    m = Model.__new__(Model)
    m.model_id, m.slug, m.verbose = model_id, "molecular-weight", False
    return m


def _defaults(fn):
    return {
        k: v.default
        for k, v in inspect.signature(fn).parameters.items()
        if v.default is not inspect.Parameter.empty
    }


def _cli_defaults(cmd):
    return {p.name: p.default for p in cmd.params}


# CLI option -> API parameter, for every option the API mirrors.
PARITY = [
    (
        serve_cmd,
        Model.serve,
        {
            "port": "port",
            "track": "track",
            "tracking_use_case": "tracking_use_case",
            "enable_cache": "enable_cache",
            "read_store": "read_store",
            "write_store": "write_store",
            "access": "access",
            "nearest_neighbors": "nearest_neighbors",
            "max_memory": "max_cache_memory_frac",
        },
    ),
    (run_cmd, Model.run, {"output": "output", "batch_size": "batch_size"}),
    (example_cmd, Model.example, {"mode": "mode", "output_file": "output"}),
    (info_cmd, Model.info, {"output": "output"}),
    (
        fetch_cmd,
        Model.fetch,
        {
            "from_dir": "from_dir",
            "from_github": "from_github",
            "from_s3": "from_s3",
            "from_hosted": "from_hosted",
            "version": "version",
        },
    ),
    (catalog_cmd, Catalog.hub, {"more": "more", "task": "task", "output": "output"}),
]


@pytest.mark.parametrize("cmd, method, mapping", PARITY)
def test_api_mirrors_cli_options_and_defaults(cmd, method, mapping):
    cli = _cli_defaults(cmd())
    api = _defaults(method)
    for cli_name, api_name in mapping.items():
        assert cli_name in cli, cli_name
        assert api_name in api, api_name
        assert api[api_name] == cli[cli_name], (cli_name, api_name)


def test_quiet_by_default_and_cli_lines_when_verbose(capsys):
    with library_call(verbose=False):
        echo("hidden", fg="green")
    assert capsys.readouterr().out == ""
    with library_call(verbose=True):
        echo("shown", fg="green")
    assert "shown" in capsys.readouterr().out


def test_decorated_errors_are_raised_not_exited(capsys):
    @throw_ersilia_exception()
    def fails():
        raise InvalidModelIdentifierError("eos9zzz")

    with library_call():
        with pytest.raises(InvalidModelIdentifierError):
            fails()
    assert capsys.readouterr().out == ""


def test_prompts_take_their_default_in_library_mode(monkeypatch):
    asked = MagicMock()
    monkeypatch.setattr(terminal, "raw_input_with_timeout", asked)
    with library_call():
        assert terminal.yes_no_input("Fetch it? [Y/n]", default_answer="n") is False
    asked.assert_not_called()


@pytest.mark.parametrize(
    "result, expected",
    [
        (FetchResult(True, "Model fetched successfully"), True),
        (FetchResult(False, "Model already exists on your system. ..."), False),
        (FetchResult(True, "Model eos3b5e is already available locally ..."), False),
    ],
)
def test_fetch_returns_a_bool(result, expected):
    fetcher = MagicMock()
    fetcher.return_value.fetch = MagicMock(return_value=_coro(result))
    with patch("ersilia.hub.fetch.fetch.ModelFetcher", fetcher):
        assert _model().fetch() is expected


def test_fetch_failure_raises():
    fetcher = MagicMock()
    fetcher.return_value.fetch = MagicMock(
        return_value=_coro(FetchResult(False, "Docker is not running."))
    )
    with patch("ersilia.hub.fetch.fetch.ModelFetcher", fetcher):
        with pytest.raises(ModelFetchError) as e:
            _model().fetch()
    assert e.value.reason == "Docker is not running."
    assert isinstance(e.value, ErsiliaError)


async def _return(value):
    return value


def _coro(value):
    return _return(value)


class _FakeServed:
    def __init__(self, rows=2):
        self.rows = rows
        self.calls = []

    def run(self, input, output, batch_size):
        self.calls.append((input, output, batch_size))
        n = len(pd.read_csv(input))
        pd.DataFrame({"key": ["k"] * n, "input": ["x"] * n, "mw": [1.0] * n}).to_csv(
            output, index=False
        )
        return output


def _served(m, fake):
    return (
        patch.object(Model, "_require_served", lambda self: None),
        patch.object(Model, "_served_model", lambda self: fake),
    )


def _api_tmp_dirs():
    return set(glob.glob(os.path.join(tempfile.gettempdir(), "ersilia-api-*")))


def test_run_accepts_a_list_or_csv_and_returns_a_dataframe_or_path(tmp_path):
    m, fake = _model(), _FakeServed()
    before = _api_tmp_dirs()
    csv = tmp_path / "in.csv"
    pd.DataFrame({"smiles": ["CCO", "CC", "C"]}).to_csv(csv, index=False)
    p1, p2 = _served(m, fake)
    with p1, p2:
        assert m.run(["CCO", "CC"]).shape == (2, 3)
        assert m.run(str(csv)).shape == (3, 3)
        out = tmp_path / "out.csv"
        assert m.run(str(csv), output=str(out)) == str(out) and out.exists()
    assert fake.calls[0][2] == 100
    assert _api_tmp_dirs() == before


@pytest.mark.parametrize(
    "kwargs",
    [
        {"input": "input.txt"},
        {"input": ["CCO"], "output": "out.txt"},
        {"input": ["CCO"], "output": "eos42ez_out.csv"},
    ],
)
def test_run_rejects_what_the_cli_rejects(kwargs):
    m = _model()
    p1, p2 = _served(m, _FakeServed())
    with p1, p2, pytest.raises(InvalidOptionError):
        m.run(**kwargs)


def test_run_needs_this_model_served():
    m = _model()
    with patch.object(Model, "_served_here", lambda self: False):
        with pytest.raises(ModelNotServedError):
            m.run(["CCO"])


def test_serve_refuses_while_another_model_is_served():
    m = _model()
    session = MagicMock()
    session.current_model_id.return_value = "eos42ez"
    with patch.object(Model, "_session", lambda self: session):
        with pytest.raises(SessionBusyError):
            m.serve()


@pytest.mark.parametrize(
    "kwargs",
    [
        {"access": "shared"},
        {"tracking_use_case": "production"},
        {"max_cache_memory_frac": 1.5},
        {"write_store": True},
    ],
)
def test_serve_validates_options_like_the_cli(kwargs):
    with pytest.raises(InvalidOptionError):
        _model().serve(**kwargs)


def test_example_uses_this_model_not_the_served_one(tmp_path):
    seen = {}

    class FakeGenerator:
        def __init__(self, model_id):
            seen["model_id"] = model_id

        def example(self, n_samples, file_name, mode):
            seen["mode"] = mode
            pd.DataFrame({"input": ["CCO"] * n_samples}).to_csv(file_name, index=False)

    with (
        patch.object(Model, "_require_fetched", lambda self: None),
        patch("ersilia.io.input.ExampleGenerator", FakeGenerator),
    ):
        assert _model("eos3b5e").example(2, mode="Curated") == ["CCO", "CCO"]
    assert seen == {"model_id": "eos3b5e", "mode": "curated"}


def test_example_rejects_unknown_modes():
    with pytest.raises(InvalidOptionError):
        _model().example(mode="predefined-ish")


def test_card_of_unknown_model_raises():
    with patch("ersilia.hub.content.card.ModelCard") as card:
        card.return_value.get.return_value = None
        with pytest.raises(InvalidModelIdentifierError):
            Catalog().card("eos9zzz")


def test_catalog_validates_task_and_output():
    with pytest.raises(InvalidOptionError):
        Catalog().hub(task="Prediction")
    with pytest.raises(InvalidOptionError):
        Catalog().local(output="models.txt")


def test_importing_the_api_changes_nothing_in_the_process():
    code = (
        "import asyncio, os\n"
        "os.umask(0o022)\n"
        "import ersilia.api\n"
        "assert os.umask(0o022) == 0o022, 'umask changed'\n"
        "assert 'CUDA_VISIBLE_DEVICES' not in os.environ, 'CUDA_VISIBLE_DEVICES set'\n"
        "assert not hasattr(asyncio, '_nest_patched'), 'asyncio patched'\n"
    )
    env = {k: v for k, v in os.environ.items() if k != "CUDA_VISIBLE_DEVICES"}
    result = subprocess.run(
        [sys.executable, "-c", code], env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr[-500:]


def test_importing_the_api_does_not_import_the_cli():
    code = "import sys, ersilia.api\nassert 'ersilia.cli' not in sys.modules\n"
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr[-500:]
