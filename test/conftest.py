import pytest


@pytest.fixture(autouse=True)
def _fresh_echo_state():
    # Each test starts like a new command: output on stdout until an error.
    from ersilia.utils.echo import reset

    reset()
    yield
