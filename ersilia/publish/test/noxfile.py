import os
import nox

@nox.session
def test_model_image(session):
    """Run ersilia test for model image in a controlled env."""
    # Install test harness (ersilia + testing extras)
    session.install("git+https://github.com/ersilia-os/ersilia.git#egg=ersilia[test]")

    model_id = os.environ.get("MODEL_ID")
    version = os.environ.get("MODEL_VERSION")

    if not model_id or not version:
        session.error("MODEL_ID and MODEL_VERSION must be set")

    # This calls  test/runner stack.
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
