import os
import asyncio
from ersilia.hub.fetch.fetch import ModelFetcher
from ersilia import ErsiliaModel

MODELS = ["eos0t01", "eos3b5e", "eos0t03", "eos0t04"]


def test_model_1():
    return
    MODEL_ID = MODELS[0]
    INPUT = ""
    ModelFetcher().fetch(MODEL_ID)
    em = ErsiliaModel(MODEL_ID)
    em.serve()
    em.invert(INPUT)
    em.shuffle(INPUT)
    em.close()
    assert 1 == 1


def test_model_2():
    MODEL_ID = MODELS[1]
    INPUT = "CCCC"
    mf = ModelFetcher(overwrite=True)
    asyncio.run(mf.fetch(
        MODEL_ID
    ))
    em = ErsiliaModel(MODEL_ID)
    em.serve()
    em.close()
    assert 1 == 1


def test_model_3():
    return
    MODEL_ID = MODELS[2]
    INPUT = "CCCNC"
    ModelFetcher().fetch(MODEL_ID)
    em = ErsiliaModel(MODEL_ID)
    em.serve()
    em.predict(INPUT)
    em.close()
    assert 1 == 1


def test_model_4():
    return
    MODEL_ID = MODELS[3]
    INPUT = ["CCCCN", "CCCOCC"]
    ModelFetcher().fetch(MODEL_ID)
    em = ErsiliaModel(MODEL_ID)
    em.serve()
    em.calculate(INPUT)
    em.close()
    assert 1 == 1
