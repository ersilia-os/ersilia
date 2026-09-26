import asyncio
import re
from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner

import ersilia.utils.terminal as terminal
from ersilia.cli.commands.close import close_cmd
from ersilia.cli.commands.delete import delete_cmd
from ersilia.cli.commands.info import info_cmd
from ersilia.cli.commands.run import run_cmd
from ersilia.cli.messages import no_model_served, wrong_extension
from ersilia.core.session import Session
from ersilia.utils.echo import echo
from ersilia.utils.exceptions_utils.exceptions import InvalidModelIdentifierError
from ersilia.utils.exceptions_utils.throw_ersilia_exception import (
    throw_ersilia_exception,
)
from ersilia.utils.logging import logger

ROOT = Path(__file__).resolve().parents[2] / "ersilia"


@pytest.fixture(autouse=True)
def default_verbosity():
    # As the CLI does by default: log records go to files, not the terminal.
    logger.set_verbosity(0)


@pytest.mark.parametrize(
    "fg, icon",
    [("green", "✓"), (None, "▪"), ("cyan", "▪"), ("yellow", "⚠"), ("red", "✖")],
)
def test_echo_icon_by_kind(capsys, fg, icon):
    echo("Hello.", fg=fg, bold=True)
    assert capsys.readouterr().out == f"  {icon}  Hello.\n"


def test_echo_continuation_lines_hang_under_the_text(capsys):
    echo("First line.\nSecond line.", fg="red")
    assert capsys.readouterr().out == "  ✖  First line.\n     Second line.\n"


def test_echo_does_not_render_emoji_shortcodes(capsys):
    echo("Done :thumbs_up:", fg="green")
    assert capsys.readouterr().out == "  ✓  Done :thumbs_up:\n"


def test_echo_can_print_to_stderr(capsys):
    echo("Bad.", fg="red", err=True)
    captured = capsys.readouterr()
    assert captured.out == "" and captured.err == "  ✖  Bad.\n"


def test_shared_wordings_are_errors_that_exit_1(capsys):
    with pytest.raises(SystemExit) as e:
        no_model_served()
    assert e.value.code == 1
    with pytest.raises(SystemExit) as e:
        wrong_extension([".json", ".csv"])
    assert e.value.code == 1
    assert capsys.readouterr().out == (
        "  ✖  No model is being served in this terminal.\n"
        "  ▪  Serve one first with 'ersilia serve MODEL'.\n"
        "  ✖  The output file must end in .json or .csv.\n"
    )
    no_model_served(fg="yellow")  # a warning (e.g. close): no exit


@pytest.mark.parametrize(
    "default, choices", [("n", "[y/N, No in 5 s]"), ("Y", "[Y/n, Yes in 5 s]")]
)
def test_prompt_shows_its_real_default_and_timeout(monkeypatch, default, choices):
    seen = {}

    def fake_input(prompt, default_answer, timeout):
        seen["prompt"] = prompt
        return ""

    monkeypatch.setattr(terminal, "raw_input_with_timeout", fake_input)
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    answer = terminal.yes_no_input("Fetch it now? [Y/n]", default_answer=default)
    assert seen["prompt"] == f"  ▪  Fetch it now? {choices} "
    assert answer is (default == "Y")


def test_prompt_timeout_says_what_was_chosen(monkeypatch, capsys):
    monkeypatch.setattr(terminal, "raw_input_with_timeout", lambda **kw: None)
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    assert terminal.yes_no_input("Fetch it now?", default_answer="n") is False
    assert "No answer; continuing with No." in capsys.readouterr().out


def test_prompts_without_a_terminal_use_the_default_at_once(monkeypatch, capsys):
    import ersilia.utils.echo as echo_module

    asked = []
    monkeypatch.setattr(
        terminal, "raw_input_with_timeout", lambda **kw: asked.append(1)
    )
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    assert terminal.yes_no_input("Fetch it now?", default_answer="n") is False
    assert echo_module.confirm("Continue?", default=True) is True
    assert asked == []
    out = capsys.readouterr().out
    assert "Fetch it now? No (no terminal to ask; using the default)." in out
    assert "Continue? Yes (no terminal to ask; using the default)." in out


def test_errors_show_message_and_hint_not_the_exception_block(capsys):
    @throw_ersilia_exception()
    def serve():
        raise InvalidModelIdentifierError("eos9zzz")

    with pytest.raises(SystemExit):
        serve()
    out = capsys.readouterr().out
    assert out.startswith(
        "  ✖  Model eos9zzz was not found in the Ersilia Model Hub.\n  ▪  "
    )
    assert "Ersilia exception class" not in out


def test_errors_raised_inside_coroutines_are_formatted(capsys):
    @throw_ersilia_exception()
    async def fetch():
        raise InvalidModelIdentifierError("eos9zzz")

    with pytest.raises(SystemExit):
        asyncio.run(fetch())
    assert "  ✖  Model eos9zzz was not found" in capsys.readouterr().out


def _invoke(cmd_factory, args, model_id=None):
    with (
        patch.object(Session, "current_model_id", return_value=model_id),
        patch.object(Session, "current_service_class", return_value=None),
        patch.object(Session, "current_output_source", return_value=None),
    ):
        return CliRunner().invoke(cmd_factory(), args)


def test_no_model_served_reads_the_same_everywhere():
    result = _invoke(run_cmd, ["-i", "in.csv", "-o", "out.csv"])
    assert result.output == (
        "  ✖  No model is being served in this terminal.\n"
        "  ▪  Serve one first with 'ersilia serve MODEL'.\n"
    )
    assert result.exit_code == 1


def test_info_without_a_model_points_to_the_card():
    result = _invoke(info_cmd, [])
    assert result.output == (
        "  ✖  No model is being served in this terminal.\n"
        "  ▪  Serve one first with 'ersilia serve MODEL'.\n"
        "  ▪  To read a model's card without serving it, use 'ersilia catalog --card MODEL'.\n"
    )
    assert result.exit_code == 1


def test_close_without_a_model_has_nothing_to_do():
    result = _invoke(close_cmd, [])
    assert result.output == (
        "  ⚠  No model is being served in this terminal, so there is nothing to close.\n"
    )
    assert result.exit_code == 0


def test_run_rejects_unsupported_output_extension(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "in.csv").write_text("smiles\nCCO\n")
    result = _invoke(run_cmd, ["-i", "in.csv", "-o", "out.txt"], model_id="eos42ez")
    assert result.output == "  ✖  The output file must end in .csv or .h5.\n"
    assert result.exit_code == 1


def test_delete_without_a_model_is_a_usage_error():
    result = CliRunner().invoke(delete_cmd(), [])
    assert result.exit_code == 2
    text = " ".join(re.sub(r"\x1b\[[0-9;]*m|[│╭╮╰╯─]", " ", result.output).split())
    assert "Give a model to delete, or use --all to delete every local model." in text


def test_user_commands_have_no_emoji_in_messages():
    # Icons come from echo; messages themselves carry no emoji or :shortcodes:.
    files = list((ROOT / "cli" / "commands").glob("*.py")) + [
        ROOT / "cli" / "messages.py",
        ROOT / "hub" / "fetch" / "fetch.py",
        ROOT / "hub" / "delete" / "delete.py",
        ROOT / "serve" / "autoservice.py",
    ]
    emoji = re.compile(r"[\U0001F300-\U0001FAFF☀-⛿✅]|:[a-z_]+:\"")
    offenders = [
        f"{f.name}:{i}"
        for f in files
        for i, line in enumerate(f.read_text().splitlines(), 1)
        if emoji.search(line)
    ]
    assert offenders == []


def _terminal_console(width=100):
    from io import StringIO

    from rich.console import Console

    return Console(
        file=StringIO(),
        force_terminal=True,
        color_system="standard",
        width=width,
        highlight=False,
    )


def _has_bold(ansi):
    return any("1" in code.split(";") for code in re.findall(r"\x1b\[([0-9;]*)m", ansi))


def test_running_spinner_sits_in_the_icon_column():
    import ersilia.utils.echo as echo_module

    c = _terminal_console()
    c.print(echo_module._SpinnerLine("Starting model eos3b5e"))
    line = re.sub(r"\x1b\[[0-9;]*m", "", c.file.getvalue()).rstrip("\n")
    assert re.fullmatch(r"  [⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏]  Starting model eos3b5e", line)
    assert not _has_bold(c.file.getvalue())


def test_confirm_prompts_look_like_other_lines(monkeypatch):
    import click

    import ersilia.utils.echo as echo_module

    seen = {}
    monkeypatch.setattr(
        click, "confirm", lambda text, **kw: seen.update(text=text, **kw) or True
    )
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    assert echo_module.confirm("Continue?") is True
    assert seen["text"] == "  ▪  Continue?" and seen["prompt_suffix"] == " "


def test_serve_panel_links_a_browsable_url_and_the_docs(monkeypatch):
    import ersilia.utils.terminal as terminal_module
    from ersilia.utils.terminal import print_serve_summary

    c = _terminal_console()
    monkeypatch.setattr(terminal_module, "console", c)
    print_serve_summary(
        "eos3b5e",
        "molecular-weight",
        "http://0.0.0.0:5000",
        -1,
        "pulled_docker",
        "/s",
        ["run"],
        "Disabled",
        False,
        False,
        None,
    )
    targets = re.findall(r"\x1b\]8;[^;]*;([^\x1b]+)\x1b\\", c.file.getvalue())
    assert targets == ["http://127.0.0.1:5000", "http://127.0.0.1:5000/docs"]


def test_options_accept_dashes_or_underscores_without_listing_both():
    from ersilia.cli.commands import ersilia_cli
    from ersilia.cli.commands.fetch import fetch_cmd

    fetch_cmd()
    both = CliRunner().invoke(
        ersilia_cli, ["fetch", "eos3b5e", "--from-github", "--from_s3"]
    )
    text = " ".join(re.sub(r"\x1b\[[0-9;]*m|[│╭╮╰╯─]", " ", both.output).split())
    assert "Choose only one source; got --from_github , --from_s3" in text
    help_text = CliRunner().invoke(ersilia_cli, ["fetch", "-h"]).output
    assert "--from_github" in help_text and "--from-github" not in help_text


def test_top_level_help_aligns_commands_with_options():
    from ersilia.cli.create_cli import create_ersilia_cli

    out = CliRunner().invoke(create_ersilia_cli(), ["--help"], terminal_width=100)
    lines = re.sub(r"\x1b\[[0-9;]*m", "", out.output).splitlines()
    option = next(line for line in lines if "--version" in line)
    command = next(line for line in lines if line.startswith("│ catalog"))
    assert option.index("Show the version") == command.index("List a catalog")
