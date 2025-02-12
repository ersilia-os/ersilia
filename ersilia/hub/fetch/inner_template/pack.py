# ruff: noqa

import os

from src.service import CHECKPOINTS_BASEDIR, FRAMEWORK_BASEDIR, Service, load_model


def main():
    """
    Main function to load the model, pack it into a service, and save the service.
    """
    root = os.path.dirname(os.path.realpath(__file__))
    mdl = load_model(
        os.path.join(root, "model", FRAMEWORK_BASEDIR),
        os.path.join(root, "model", CHECKPOINTS_BASEDIR),
    )

    service = Service()
    service.pack("model", mdl)
    service.save()


if __name__ == "__main__":
    main()
