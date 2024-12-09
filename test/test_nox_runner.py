import os
import pytest
from .playground.runner import NoxAPI


@pytest.fixture
def config_path():
    """Fixture to provide the path to the Nox configuration file."""
    return os.path.join(os.getcwd(), "playground/config.yml")


def test_execute_sessions(config_path):
    """Test executing Nox sessions."""
    # Initialize the NoxAPI with the given configuration path
    api = NoxAPI(config_path=config_path)

    # Add sessions to the queue
    api.setup()
    api.test_from_github()

    # Execute all sessions
    api.execute_all()

    # Validate the queue is cleared after execution
    assert len(api.queue) == 0

    # Check that expected session outputs exist (logs or side effects)
    # Example: Verify that `setup` and `test_from_github` ran successfully
    # # Replace this with project-specific checks
    # session_log_path = os.path.join(os.getcwd(), "nox_output.log")
    # with open(session_log_path, "r") as log_file:
    #     log_content = log_file.read()
    #     assert "Setup session executed." in log_content
    #     assert "Test from GitHub session executed." in log_content
