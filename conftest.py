import io
from contextlib import redirect_stdout
from test.playground.shared import results


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    command_col_width = 50
    description_col_width = 15
    time_col_width = 15
    memory_col_width = 15
    status_col_width = 20

    border = (
        f"+{'-' * command_col_width}+{'-' * description_col_width}+"
        f"{'-' * time_col_width}+{'-' * memory_col_width}+{'-' * status_col_width}+"
    )

    terminalreporter.write("\n\nCommand Execution Summary:\n")
    terminalreporter.write(f"{border}\n")
    terminalreporter.write(
        f"| {'Command':^{command_col_width}} | {'Description':^{description_col_width}} | "
        f"{'Time Taken':^{time_col_width}} | {'Max Memory':^{memory_col_width}} | {'Status':^{status_col_width}} \n"
    )
    terminalreporter.write(f"{border}\n")

    for result in results:
        terminalreporter.write(
            f"| {result['command']:<{command_col_width}} | {result['description']:<{description_col_width}} | "
            f"{result['time_taken']:<{time_col_width}} | {result['max_memory']:<{memory_col_width}} | "
            f"{result['status']:<{status_col_width}} |\n"
        )

    terminalreporter.write(f"{border}\n")
    terminalreporter.write("\n")
