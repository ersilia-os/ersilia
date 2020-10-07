"""
from bentoml.cli import create_bentoml_cli

def create_ersilia_cli():
    _cli = create_bentoml_cli()
    # TODO Check BentoML example to produce something similar.

cli = create_ersilia_cli()

if __name__ == "__main__":
    cli()
"""
from .ersilia_service import create_ersilia_service_cli

cli = create_ersilia_service_cli()

if __name__ == "__main__":
    cli()