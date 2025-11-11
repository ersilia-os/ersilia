import os
import nox

@nox.session(venv_backend='conda', python='3.12')
def test_model_image(session):
    session.env["ERSILIA_NOX"] = "1"
#    session.env["CONDA_NO_PLUGINS"] = "true"   
    session.install("git+https://github.com/ersilia-os/ersilia.git#egg=ersilia[test]")
#    session.install('-e',"../../../.[test]") #For local testing
    model_id = os.environ.get("MODEL_ID")
    version = os.environ.get("MODEL_VERSION")
    if not model_id or not version:
        session.error("MODEL_ID and MODEL_VERSION must be set")

    session.run(
        "ersilia",
        "-v",
        "test",
        model_id,
        "--deep",
        "--from_dockerhub",
        "--version",
        version,
        external=True,
    )
