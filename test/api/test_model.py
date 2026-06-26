from ersilia.api import Model

MODEL_ID = "eos3b5e"

INPUT_LIST = [
    "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O"
]

def test_fetch():
    mdl = Model(MODEL_ID)
    assert mdl.fetch()
    assert mdl.is_fetched()

def test_serve_run_and_close():
    mdl = Model(MODEL_ID)
    serve_result = mdl.serve()
    assert serve_result

    df = mdl.run(input_list = INPUT_LIST)
    assert df.shape[0] == len(INPUT_LIST)

    close_result = mdl.close()
    assert close_result

def test_delete():
    mdl = Model(MODEL_ID)
    assert mdl.delete()
