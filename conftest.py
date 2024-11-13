import io
from contextlib import redirect_stdout
from test.playground.shared import results

# Helper function to wrap text to a specified width
def wrap_text(text, width):
    wrapped_lines = []
    words = text.split()
    current_line = ""

    for word in words:
        if len(current_line) + len(word) + 1 > width:
            wrapped_lines.append(current_line)
            current_line = word
        else:
            if current_line:
                current_line += " "
            current_line += word
    wrapped_lines.append(current_line)

    return wrapped_lines

def pytest_terminal_summary(terminalreporter, exitstatus, config):
    command_col_width = 50
    description_col_width = 15
    time_col_width = 15
    memory_col_width = 15
    status_col_width = 20
    checkups_col_width = 50

    border = (
        f"+{'-' * command_col_width}+{'-' * description_col_width}+"
        f"{'-' * time_col_width}+{'-' * memory_col_width}+{'-' * status_col_width}+{'-' * checkups_col_width}+"
    )

    terminalreporter.write("\n\nCommand Execution Summary:\n")
    terminalreporter.write(f"{border}\n")
    terminalreporter.write(
        f"| {'Commands':^{command_col_width}} | {'Descrpt.':^{description_col_width}} | "
        f"{'Time':^{time_col_width}} | {'Max Mem':^{memory_col_width}} | {'Status':^{status_col_width}} | "
        f"{'Checkups':^{checkups_col_width}} |\n"
    )
    terminalreporter.write(f"{border}\n")

    for result in results:
        checkups = ", ".join([f"{check['name']}: {'✔' if check['status'] else '✘'}" for check in result["checkups"]])
        
        # Wrap each column's text to fit within its width
        command_lines = wrap_text(result["command"], command_col_width)
        description_lines = wrap_text(result["description"], description_col_width)
        time_taken_lines = wrap_text(result["time_taken"], time_col_width)
        max_memory_lines = wrap_text(result["max_memory"], memory_col_width)
        status_lines = wrap_text(result["status"], status_col_width)
        checkups_lines = wrap_text(checkups, checkups_col_width)

        # Get the maximum number of lines needed to display all wrapped text for this row
        max_lines = max(
            len(command_lines), len(description_lines), len(time_taken_lines),
            len(max_memory_lines), len(status_lines), len(checkups_lines)
        )

        # Print each line for the current result, ensuring alignment across columns
        for i in range(max_lines):
            terminalreporter.write(
                f"| {command_lines[i] if i < len(command_lines) else '':<{command_col_width}} "
                f"| {description_lines[i] if i < len(description_lines) else '':<{description_col_width}} "
                f"| {time_taken_lines[i] if i < len(time_taken_lines) else '':<{time_col_width}} "
                f"| {max_memory_lines[i] if i < len(max_memory_lines) else '':<{memory_col_width}} "
                f"| {status_lines[i] if i < len(status_lines) else '':<{status_col_width}} "
                f"| {checkups_lines[i] if i < len(checkups_lines) else '':<{checkups_col_width}} |\n"
            )

        # Print a dividing line after each result entry
        terminalreporter.write(f"{border}\n")

    terminalreporter.write("\n")
