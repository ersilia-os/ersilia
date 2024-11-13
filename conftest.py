from test.playground.shared import results
from rich.table import Table
from rich.console import Console
from rich.text import Text
from rich import box

def pytest_terminal_summary(terminalreporter, exitstatus, config):
    console = Console()
    table = Table(title="Command Execution Summary", box=box.SQUARE)

    table.add_column("Command", width=50)
    table.add_column("Description", width=15)
    table.add_column("Time Taken", width=15)
    table.add_column("Max Memory", width=15)
    table.add_column("Status", width=20)
    table.add_column("Checkups", width=30)
    table.add_column("Docker Status", width=20)
    table.add_column("Runner", width=20)
    table.add_column("CLI Type", width=20)

    for result in results:
        formatted_checkups = []
        for check in result["checkups"]:
            if check["status"]:
                formatted_checkups.append(Text("✔", style="green") + f" {check['name']}")
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
            Text("✔", style="green") if result["activate_docker"] else Text("✘", style="red"),
            result["runner"],
            result["cli_type"],
        )

    console.print(table)
