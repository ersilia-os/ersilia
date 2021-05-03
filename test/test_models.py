from ersilia import ErsiliaModel


def test_eos0t01():
    MODEL_ID = "eos0t01"
    INPUT = ""
    em = ErsiliaModel(MODEL_ID)
    em.serve()
    em.invert(INPUT)
    em.shuffle(INPUT)
    em.close()
    assert 1 == 1


def test_eos0t02():
    MODEL_ID = "eos0t02"
    INPUT = "CCCC"
    em = ErsiliaModel(MODEL_ID)
    em.serve()
    em.predict(INPUT)
    em.close()
    assert 1 == 1


def test_eos0t03():
    MODEL_ID = "eos0t03"
    INPUT = "CCCNC"
    em = ErsiliaModel(MODEL_ID)
    em.serve()
    em.predict(INPUT)
    em.close()
    assert 1 == 1
