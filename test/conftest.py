import pytest


@pytest.fixture(autouse=True)
def _fresh_command_state(monkeypatch):
    # Each test starts like a new command: output on stdout until an error,
    # and no models stopped by an earlier cleanup of orphaned sessions.
    import ersilia.utils.session as session_utils
    from ersilia.utils.echo import reset

    reset()
    monkeypatch.setattr(session_utils, "_stopped_by_cleanup", [])
    yield
