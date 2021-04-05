MODEL_ID = "eos0t01"
INPUT = {}


def test_import():
    import ersilia
    assert (1 == 1)


def test_fetch_conda():
    from ersilia.hub.fetch import ModelFetcher
    mf = ModelFetcher(overwrite=False)
    mf.fetch(model_id=MODEL_ID, dockerize=False)
    assert (1 == 1)


def test_systembundle_predict():
    from ersilia.serve.services import SystemBundleService
    srv = SystemBundleService(MODEL_ID)
    srv.serve()
    output = srv.api("predict", INPUT)
    srv.close()
    assert (output == "Hello Ersilia!")
