import pytest
from ersilia.db.hubdata.interfaces import JsonModelsInterface

ji = JsonModelsInterface(config_json=None)

import pytest


class OriginalModelFinder:
    def __init__(self):
        pass

    def _read_json_file(self):
        return ji._read_json_file()

    def items(self):
        models = self._read_json_file()
        for mdl in models:
            yield mdl

    def _find_url_using_s3_models_json(self, model_id):
        for mdl in self.items():
            if mdl["Identifier"] == model_id:
                if "Host URL" in mdl:
                    return mdl["Host URL"]
                else:
                    return None
        return None


class OptimizedModelFinder:
    def __init__(self):
        self.models_cache = None

    def _read_json_file(self):
        return ji._read_json_file()

    def _cache_models(self):
        if self.models_cache is None:
            models = self._read_json_file()
            self.models_cache = {mdl["Identifier"]: mdl for mdl in models}

    def _find_url_using_s3_models_json(self, model_id):
        self._cache_models()
        model = self.models_cache.get(model_id)

        if model:
            if "Host URL" in model:
                return model["Host URL"]
            else:
                return None
        return None


@pytest.fixture
def original_finder():
    return OriginalModelFinder()


@pytest.fixture
def actual_url():
    return "https://eos9ei3-tkreo.ondigitalocean.app/"


@pytest.fixture
def optimized_finder():
    return OptimizedModelFinder()


def test_original_finder(actual_url, original_finder):
    url = original_finder._find_url_using_s3_models_json("eos9ei3")
    assert url == actual_url


def test_optimized_finder(actual_url, optimized_finder):
    url = optimized_finder._find_url_using_s3_models_json("eos9ei3")
    assert url == actual_url


def test_benchmark_original(actual_url, benchmark, original_finder):
    result = benchmark(original_finder._find_url_using_s3_models_json, "eos9ei3")
    assert result == actual_url


def test_benchmark_optimized(actual_url, benchmark, optimized_finder):
    result = benchmark(optimized_finder._find_url_using_s3_models_json, "eos9ei3")
    assert result == actual_url
