from .bentoml import BentoMLRequirement


def check_bentoml():
    req = BentoMLRequirement()
    if not req.is_bentoml_ersilia_version():
        req.install()
    return
