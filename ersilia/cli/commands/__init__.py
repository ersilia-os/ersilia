from pathlib import Path

import rich_click as click
from rich_click import RichCommand, RichGroup

from ... import __version__ as __version__
from ... import logger

_ASCII_ART_FILE = Path(__file__).parents[3] / "assets" / "ascii-art.txt"


def _print_logo():
    if not _ASCII_ART_FILE.exists():
        return
    from rich.console import Console
    from rich.text import Text

    # Ersilia palette — deeper tones visible on both dark and light terminals
    gradient = ["#E07040", "#C060C0", "#8855BB", "#5599DD", "#55BB77"]

    console = Console(highlight=False)
    lines = _ASCII_ART_FILE.read_text().splitlines()
    non_empty = [l for l in lines if l.strip()]
    n = max(len(non_empty) - 1, 1)
    ci = 0
    for line in lines:
        if line.strip():
            color = gradient[ci * (len(gradient) - 1) // n]
            console.print(Text(line, style=f"bold {color}"), end="\n")
            ci += 1
        else:
            console.print("")


click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.COLOR_SYSTEM = "truecolor"
click.rich_click.STYLE_OPTION = "bold magenta"
click.rich_click.STYLE_COMMAND = "bold green"
click.rich_click.STYLE_METAVAR = "italic yellow"
click.rich_click.STYLE_SWITCH = "underline cyan"
click.rich_click.STYLE_USAGE = "bold blue"
click.rich_click.STYLE_OPTION_DEFAULT = "dim italic"
click.rich_click.STYLE_HELPTEXT = ""
click.rich_click.STYLE_ERRORS_SUGGESTION = "bold"
click.rich_click.HEADER_TEXT = f"Ersilia version: {__version__}"


def _align_commands_with_options():
    # rich-click stretches the commands table to the panel width and gives
    # the extra space to the name column, so command descriptions start
    # further right than option descriptions. Size the name column to the
    # option columns instead (options here are flags, with no metavar).
    try:
        from rich_click.rich_panel import RichCommandPanel
    except ImportError:
        return
    get_table = RichCommandPanel.get_table

    def aligned_get_table(self, command, ctx, formatter):
        table = get_table(self, command, ctx, formatter)
        try:
            names = [
                p.opts + p.secondary_opts
                for p in command.get_params(ctx)
                if isinstance(p, click.Option) and not p.hidden
            ]
            long_ = max(
                len("/".join(o for o in n if o.startswith("--"))) for n in names
            )
            short = max(
                len("/".join(o for o in n if not o.startswith("--"))) for n in names
            )
            table.expand = False
            table.columns[0].min_width = long_ + short + 2 if short else long_ + 1
        except Exception:
            pass
        return table

    RichCommandPanel.get_table = aligned_get_table


_align_commands_with_options()

# ruff: noqa: D101, D102


class ErsiliaCommandGroup(RichGroup):
    def command(self, *args, **kwargs):
        kwargs.setdefault("cls", RichCommand)
        return RichGroup.command(self, *args, **kwargs)

    def main(self, *args, **kwargs):
        import sys

        argv = sys.argv[1:]
        if not argv or "--help" in argv or "-h" in argv:
            _print_logo()
        return super().main(*args, **kwargs)

    def invoke(self, ctx):
        # Ctrl+C anywhere (including at a prompt) ends with one short line,
        # never a traceback.
        try:
            return super().invoke(ctx)
        except (KeyboardInterrupt, click.exceptions.Abort) as e:
            import sys

            from ..echo import echo

            echo("Interrupted.", fg="yellow")
            note = getattr(e, "ersilia_note", None) or getattr(
                e.__context__, "ersilia_note", None
            )
            if note:
                echo(note)
            sys.exit(130)

    def resolve_command(self, ctx, args):
        # Suggest the closest command for a typo, e.g. 'ersilia server'.
        name = args[0] if args else None
        if name and not name.startswith("-") and self.get_command(ctx, name) is None:
            import difflib

            matches = difflib.get_close_matches(name, self.list_commands(ctx), n=1)
            if matches:
                ctx.fail(f"No such command '{name}'. Did you mean '{matches[0]}'?")
        return super().resolve_command(ctx, args)


@click.group(
    cls=ErsiliaCommandGroup,
    context_settings={
        "show_default": True,
        "help_option_names": ["-h", "--help"],
        # --from-github works like --from_github, without listing both
        # spellings in the help.
        "token_normalize_func": lambda name: name.replace("-", "_"),
    },
    epilog="To learn more about a specific command, run: ersilia COMMAND --help",
)
@click.version_option(version=__version__)
@click.option(
    "-v",
    "--verbose",
    default=False,
    is_flag=True,
    help="Show logging on terminal when running commands.",
)
def ersilia_cli(verbose):
    if verbose:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)
