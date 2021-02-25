from ersilia import ErsiliaModel

inputs = [
    "CCCC"
]

def predict_from_pip():
    em = ErsiliaModel("eos0aaa")
    em.predict(inputs)


def predict_from_conda():
    em = ErsiliaModel("eos0aaa")
    em.predict(inputs)


def predict_from_docker():
    em = ErsiliaModel("eos0aaa")
    em.predict(inputs)


def predict_from_remote():
    em = ErsiliaModel("eos0aaa")
    em.predict(inputs)
