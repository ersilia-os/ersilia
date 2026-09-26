"""Integration test of the Python API catalog: needs network access."""

from ersilia.api import Catalog


def test_catalog_hub():
    df = Catalog().hub()
    assert df.shape[0] > 10
    assert "Identifier" in df.columns


def test_catalog_hub_by_task():
    everything = Catalog().hub()
    sampling = Catalog().hub(task="Sampling")
    assert 0 < sampling.shape[0] < everything.shape[0]


def test_catalog_local_is_a_dataframe():
    df = Catalog().local()
    assert list(df.columns)


def test_catalog_card():
    card = Catalog().card("eos3b5e")
    assert card["Identifier"] == "eos3b5e"
