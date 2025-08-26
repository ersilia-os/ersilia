from ersilia.api import Model

MODEL_ID = "eos3b5e"

INPUT_LIST = [
    "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O"
]

def test_fetch():
    mdl = Model(MODEL_ID)
    mdl.fetch()
    assert mdl.is_fetched()

def test_serve_run_and_close():
    mdl = Model(MODEL_ID)
    mdl.serve()
    df = mdl.run(input_list = INPUT_LIST)
    assert df.shape[0] == len(INPUT_LIST)
    mdl.close()
    assert 1 == 1

def test_delete():
    mdl = Model(MODEL_ID)
    mdl.delete()
    assert 1 == 1
