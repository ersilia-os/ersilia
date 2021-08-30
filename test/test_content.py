from ersilia.hub.content.catalog import ModelCatalog
from ersilia.hub.content.card import ModelCard

MODEL_ID = "eos0t01"


def test_catalog_local():
    mc = ModelCatalog()
    mc.local()
    assert 1 == 1


def test_catalog_hub():
    mc = ModelCatalog()
    mc.hub()
    assert 1 == 1


def test_card():
    mc = ModelCard()
    mc.get(MODEL_ID)
    assert 1 == 1
