"""Pull progress from Docker Engine API events; no Docker needed."""

from unittest.mock import MagicMock, patch

import pytest

from ersilia.hub.pull.pull import PullProgress, pull_with_progress

MB = 1_000_000


def _ev(layer, status, current=None, total=None):
    event = {"id": layer, "status": status}
    if total is not None:
        event["progressDetail"] = {"current": current, "total": total}
    return event


def _feed(progress, events):
    for event in events:
        progress.update(event)
    return progress


def test_reports_megabytes_and_layers_while_downloading():
    p = _feed(
        PullProgress(expected_bytes=30 * MB),
        [
            {"id": "latest", "status": "Pulling from ersiliaos/eos3b5e"},
            _ev("aaa", "Pulling fs layer"),
            _ev("bbb", "Pulling fs layer"),
            _ev("aaa", "Downloading", 5 * MB, 10 * MB),
            _ev("bbb", "Downloading", 10 * MB, 20 * MB),
        ],
    )
    assert p.describe() == "15/30 MB  0/2 layers"
    completed, total = p.bar()
    assert (completed, total) == (15 * MB, 60 * MB)


def test_switches_to_extracting_once_everything_is_downloaded():
    p = _feed(
        PullProgress(),
        [
            _ev("aaa", "Downloading", 5 * MB, 10 * MB),
            _ev("bbb", "Downloading", 20 * MB, 20 * MB),
            _ev("bbb", "Extracting", 5 * MB, 20 * MB),
        ],
    )
    assert not p.extracting  # aaa is still downloading
    _feed(
        p, [_ev("aaa", "Download complete"), _ev("bbb", "Extracting", 10 * MB, 20 * MB)]
    )
    assert p.describe() == "extracting 10/30 MB  0/2 layers"
    _feed(p, [_ev("bbb", "Pull complete"), _ev("aaa", "Pull complete")])
    assert p.describe() == "30/30 MB  2/2 layers"
    assert p.bar() == (60 * MB, 60 * MB)


def test_layers_already_present_are_not_counted_as_download():
    p = _feed(
        PullProgress(expected_bytes=100 * MB),
        [
            _ev("aaa", "Already exists"),
            _ev("bbb", "Downloading", 1 * MB, 4 * MB),
        ],
    )
    assert p.describe() == "1/4 MB  1/2 layers"


def test_without_byte_counts_the_bar_follows_layers():
    # Some daemons report statuses only; the bar must still move.
    p = _feed(
        PullProgress(),
        [
            {"id": "aaa", "status": "Downloading"},
            {"id": "bbb", "status": "Downloading"},
            {"id": "aaa", "status": "Pull complete"},
        ],
    )
    assert p.bar() == (1, 2)
    assert p.describe() == "1/2 layers"


def test_nothing_known_yet_means_an_indeterminate_bar():
    assert PullProgress().bar() == (0, None)


def test_pull_errors_are_raised():
    client = MagicMock()
    client.api.pull.return_value = iter(
        [_ev("aaa", "Downloading", 1, 2), {"error": "no matching manifest"}]
    )
    with (
        patch("docker.from_env", return_value=client),
        patch("ersilia.utils.docker.set_docker_host"),
    ):
        with pytest.raises(RuntimeError, match="no matching manifest"):
            pull_with_progress("ersiliaos/eos3b5e", "latest", lambda p: None)


def test_known_layer_sizes_give_a_fixed_total_from_the_start():
    sizes = {"aaa": 10 * MB, "bbb": 20 * MB, "ccc": 5 * MB}
    p = _feed(PullProgress(layer_sizes=sizes), [_ev("aaa", "Pulling fs layer")])
    assert p.total_bytes == 35 * MB  # bbb and ccc not announced yet, still counted
    _feed(p, [_ev("ccc", "Already exists"), _ev("bbb", "Pulling fs layer")])
    assert p.total_bytes == 30 * MB  # ccc is present locally
    _feed(p, [_ev("aaa", "Downloading", 4 * MB, 10 * MB)])
    assert p.describe() == "4/30 MB  1/3 layers"
