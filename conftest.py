from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from test.playground.shared import results


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    console = Console()

    docker_status = (
        Text("✔", style="green")
        if any(result["activate_docker"] for result in results)
        else Text("✘", style="red")
    )
    runner = results[0]["runner"] if results else "N/A"
    cli_type = results[0]["cli_type"] if results else "N/A"

    header_panel = Panel.fit(
        f"Docker Status: {docker_status}\nRunner: {runner}\nCLI Type: {cli_type}",
        title="Execution Summary",
        border_style="bold",
    )
    console.print(header_panel)

    table = Table(title="Command Execution Summary", box=box.SQUARE)
    table.add_column("Command", width=50)
    table.add_column("Description", width=15)
    table.add_column("Time Taken", width=15)
    table.add_column("Max Memory", width=15)
    table.add_column("Status", width=20)
    table.add_column("Checkups", width=30)

    for result in results:
        formatted_checkups = []
        for check in result["checkups"]:
            if check["status"]:
                formatted_checkups.append(
                    Text("✔", style="green") + f" {check['name']}"
                )
            else:
                formatted_checkups.append(Text("✘", style="red") + f" {check['name']}")
        checkups_text = "\n".join(str(checkup) for checkup in formatted_checkups)

        table.add_row(
            result["command"],
            result["description"],
            result["time_taken"],
            result["max_memory"],
            result["status"],
            checkups_text,
        )

    console.print(table)
