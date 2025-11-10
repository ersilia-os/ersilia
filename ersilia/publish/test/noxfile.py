import os
import nox

@nox.session(
  venv_backend="conda", python=["3.8", "3.9", "3.10", "3.11", "3.12"], reuse_venv=True
)
def test_model_image(session):
    session.env["ERSILIA_NOX"] = "1"  
    session.install("git+https://github.com/ersilia-os/ersilia.git#egg=ersilia[test]")

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
