import json
from ersilia.utils.terminal import truncate_output


def test_truncate_output_list():
    # Test truncation for a long list
    output = list(range(20))  # List with 20 items
    truncated = truncate_output(output, max_items=10)
    assert truncated == "[0, 1, 2, 3, 4, 5, 6, 7, 8, 9] ... (and 10 more items)"


def test_truncate_output_dict():
    # Test truncation for a long dictionary
    output = {f"key{i}": i for i in range(20)}  # Dictionary with 20 key-value pairs
    truncated = truncate_output(output, max_items=5)
    assert truncated.startswith('{\n    "key0": 0,')
    assert truncated.endswith("... (and 15 more lines)")


def test_truncate_output_short_list():
    # Test for a short list that doesn't need truncation
    output = [1, 2, 3]
    truncated = truncate_output(output, max_items=10)
    assert truncated == "[1, 2, 3]"


def test_truncate_output_short_dict():
    # Test for a short dictionary that doesn't need truncation
    output = {"key1": 1, "key2": 2}
    truncated = truncate_output(output, max_items=10)
    assert "key1" in truncated
    assert "key2" in truncated


def test_truncate_output_short_string():
    # Test for a short string that doesn't need truncation
    output = "Short string"
    truncated = truncate_output(output, max_chars=50)
    assert truncated == "Short string"


def test_truncate_output_other_types():
    # Test for non-list, non-dict, non-string types
    output = 12345
    truncated = truncate_output(output)
    assert truncated == "12345"
