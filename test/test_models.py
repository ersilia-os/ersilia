from ersilia.hub.fetch.fetch import ModelFetcher
from ersilia import ErsiliaModel

MODELS = ["eos0t01", "eos0t02", "eos0t03"]


def test_model_0():
    MODEL_ID = MODELS[0]
    INPUT = ""
    mf = ModelFetcher().fetch(MODEL_ID)
    em = ErsiliaModel(MODEL_ID)
    em.serve()
    em.invert(INPUT)
    em.shuffle(INPUT)
    em.close()
    assert 1 == 1


def test_model_1():
    MODEL_ID = MODELS[1]
    INPUT = "CCCC"
    mf = ModelFetcher().fetch(MODEL_ID)
    em = ErsiliaModel(MODEL_ID)
    em.serve()
    pred = em.predict(INPUT)
    em.close()
    assert 1 == 1


def test_model_2():
    MODEL_ID = MODELS[2]
    INPUT = "CCCNC"
    mf = ModelFetcher().fetch(MODEL_ID)
    em = ErsiliaModel(MODEL_ID)
    em.serve()
    em.predict(INPUT)
    em.close()
    assert 1 == 1
