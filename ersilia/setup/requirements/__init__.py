from .bentoml_requirement import BentoMLRequirement


def check_bentoml():
    req = BentoMLRequirement()
    if not req.is_installed():
        req.install()
    if not req.is_bentoml_ersilia_version():
        req.install()
    return
