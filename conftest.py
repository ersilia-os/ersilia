from pathlib import Path

from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from ersilia.default import EOS_PLAYGROUND
from test.playground.shared import results

base_path = Path(EOS_PLAYGROUND)
log_path = base_path / "logs"
fil_path = base_path / "files"


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    console = Console()

    if not results:
        console.print(
            Panel.fit("No results available.", title="Summary", border_style="bold")
        )
        return

    docker_status = (
        Text("✔", style="green bold")
        if any(result["activate_docker"] for result in results)
        else Text("✘", style="red bold")
    )

    cli_value = results[0]["cli"]
    if isinstance(cli_value, list):
        cli_value = ", ".join(map(str, cli_value))

    runner = Text(results[0]["runner"], style="bold") if results else "N/A"
    cli_type = Text(cli_value, style="bold") if results else "N/A"
    show_remark = bool(results[0]["show_remark"])
    header_panel = Panel.fit(
        f"Docker Status: {docker_status}\nRunner: {runner}\
            \nCli Types: {cli_type}\nLogs Path: {log_path}\nInput/Output Path: {fil_path}",
        title="Summary",
        border_style="bold",
    )
    console.print(header_panel)

    table = Table(title="Command Execution Summary", box=box.SQUARE)
    table.add_column("Command", width=40)
    table.add_column("Time Taken", width=15, justify="center")
    table.add_column("Max Memory", width=15, justify="center")
    table.add_column("Status", width=10, justify="right")
    table.add_column("Checkups", width=50, justify="left")
    if show_remark:
        table.add_column("Remark", width=70, justify="right")

    for result in results:
        formatted_checkups = []
        for check in result["checkups"]:
            if check["status"]:
                formatted_checkups.append(Text(f"✔ {check['name']}", style="green"))
            else:
                formatted_checkups.append(Text(f"✘ {check['name']}", style="red"))
        checkups_text = "\n".join(str(checkup) for checkup in formatted_checkups)

        table.add_row(
            Text(result["command"], style="bold"),
            result["time_taken"],
            result["max_memory"],
            result["status"],
            checkups_text,
            Text(result["remark"], style="bold") if show_remark else "",
        )

    console.print(table)


def pytest_collection_modifyitems(config, items):
    desired = [
        "test/cli/test_catalog.py",
        "test/cli/test_close.py",
        "test/cli/test_fetch.py",
        "test/cli/test_serve.py",
        "test/cli/test_run.py",
        "test/cli/test_delete.py",
    ]
    desired = [p.replace("\\", "/") for p in desired]
    order = {p: i for i, p in enumerate(desired)}
    original_index = {item: idx for idx, item in enumerate(items)}

    def key(it):
        file_path = it.nodeid.split("::", 1)[0].replace("\\", "/")
        return (order.get(file_path, len(desired)), original_index[it])

    items.sort(key=key)
