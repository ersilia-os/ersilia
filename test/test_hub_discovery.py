MODEL_ID = "eos0t01"


def test_catalog():
    from ersilia.hub.catalog import ModelCatalog

    mc = ModelCatalog()
    mc.local()
    assert 1 == 1


def test_card():
    from ersilia.hub.card import ModelCard

    mc = ModelCard()
    mc.get(MODEL_ID)
    assert 1 == 1
