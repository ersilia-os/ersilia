
def test_eos0t01():
    MODEL_ID = "eos0t01"
    INPUT = {}
    from ersilia import ErsiliaModel
    em = ErsiliaModel(MODEL_ID)
    em.serve()
    em.invert(INPUT)
    em.shuffle(INPUT)
    em.close()
    assert (1 == 1)
