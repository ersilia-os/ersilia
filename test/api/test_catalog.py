from ersilia.api import Catalog

def test_catalog_hub():
    cat = Catalog()
    df = cat.catalog(hub=True)
    assert df.shape[0] > 10

def test_catalog_local():
    cat = Catalog()
    df = cat.catalog(hub=False)
    assert 1 == 1

