import subprocess
import sys
import time
from unittest.mock import MagicMock

import psutil

from ersilia.serve.services import _FastApiService

# Stand-in for ersilia_model_serve: it starts the model server (run_uvicorn.py)
# as a child with subprocess.run and waits on it.
PARENT_SCRIPT = (
    "import subprocess, sys; "
    "subprocess.run([sys.executable, '-c', 'import time; time.sleep(60)'])"
)


def _service_for(pid):
    service = _FastApiService.__new__(_FastApiService)
    service.pid = pid
    service.logger = MagicMock()
    return service


def test_close_stops_server_child_process():
    """close() must stop the model server child, not only the parent (#1921)."""
    parent = subprocess.Popen([sys.executable, "-c", PARENT_SCRIPT])
    try:
        for _ in range(50):
            children = psutil.Process(parent.pid).children(recursive=True)
            if children:
                break
            time.sleep(0.1)
        assert children, "server child process did not start"

        _service_for(parent.pid).close()
        parent.wait(timeout=10)

        _, alive = psutil.wait_procs(children, timeout=10)
        assert not alive, "server child process is still running after close()"
    finally:
        for proc in psutil.Process().children(recursive=True):
            proc.kill()


def test_close_with_missing_process_does_not_raise():
    """close() on a PID that no longer exists only logs."""
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    proc.wait()
    service = _service_for(proc.pid)
    service.close()
    service.logger.info.assert_called_once()


def test_close_with_unset_pid_does_not_touch_current_process():
    """psutil.Process(None) is the current process; close() must not use it."""
    for pid in (None, 0, -1):
        service = _service_for(pid)
        service.close()
        service.logger.info.assert_called_once()
    service = _FastApiService.__new__(_FastApiService)
    service.logger = MagicMock()
    service.close()
    service.logger.info.assert_called_once()
