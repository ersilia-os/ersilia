MODEL_ID = "eos0aaa"
INPUT = ["CCC", "CCN"]

def test_import():
    import ersilia
    assert (1 == 1)

def test_fetch_conda():
    from ersilia.hub.fetch import ModelFetcher
    mf = ModelFetcher(overwrite=False)
    mf.fetch(model_id=MODEL_ID, containerize=False)
    assert (1 == 1)

def test_conda_predict():
    from ersilia.serve.services import CondaEnvironmentService
    srv = CondaEnvironmentService(MODEL_ID)
    srv.serve()
    srv.api("predict", INPUT)
    srv.close()
    assert (1 == 1)
