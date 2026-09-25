import json
import logging
import socket
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest
import requests

import ersilia.serve.standard_api as standard_api
from ersilia.serve.standard_api import StandardCSVRunApi


class _Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    mode = "ok"
    requests_seen = []
    ports_seen = set()

    def do_POST(self):
        body = self.rfile.read(int(self.headers["Content-Length"]))
        batch = json.loads(body)
        type(self).requests_seen.append(len(batch))
        type(self).ports_seen.add(self.client_address[1])
        if self.mode == "slow":
            time.sleep(2)
        if self.mode == "error":
            self.send_response(500)
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        out = json.dumps([{"value": 1} for _ in batch]).encode()
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(out)))
        self.end_headers()
        self.wfile.write(out)

    def log_message(self, *args):
        pass


@pytest.fixture
def server():
    _Handler.mode = "ok"
    _Handler.requests_seen = []
    _Handler.ports_seen = set()
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_address[1]}/run"
    httpd.shutdown()
    httpd.server_close()


def _api():
    api = StandardCSVRunApi.__new__(StandardCSVRunApi)
    api.http = requests.Session()
    api.local_cache = False
    api.logger = logging.getLogger("test_run_http")
    return api


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def test_batches_reuse_one_connection(server):
    api = _api()
    for _ in range(5):
        values, _ = api._post_batch_with_fallback(server, ["C", "CC"])
        assert values == [{"value": 1}, {"value": 1}]
    assert len(_Handler.ports_seen) == 1


def test_dead_server_fails_fast_without_splitting():
    api = _api()
    url = f"http://127.0.0.1:{_free_port()}/run"
    st = time.perf_counter()
    with pytest.raises(RuntimeError, match="did not respond"):
        api._post_batch_with_fallback(url, ["C"] * 100)
    assert time.perf_counter() - st < 5


def test_hung_server_fails_fast_without_splitting(server, monkeypatch):
    monkeypatch.setattr(standard_api, "RUN_READ_TIMEOUT", 0.5)
    _Handler.mode = "slow"
    api = _api()
    with pytest.raises(RuntimeError, match="did not respond"):
        api._post_batch_with_fallback(server, ["C"] * 8)
    assert _Handler.requests_seen == [8]


def test_server_errors_still_split_batch(server):
    _Handler.mode = "error"
    api = _api()
    values, _ = api._post_batch_with_fallback(server, ["C"] * 4)
    assert values == [None] * 4
    assert _Handler.requests_seen == [4, 2, 1, 1, 2, 1, 1]
