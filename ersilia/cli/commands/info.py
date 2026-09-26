import csv
import json

import rich_click as click

from .. import echo
from . import ersilia_cli


def info_cmd():
    """
    Provides information about a specified model.

    This command allows users to get detailed information about a current active session.

    Returns
    -------
    function
        The info command function to be used by the CLI.
    """

    @ersilia_cli.command(
        short_help="Get model information",
        help="Display information about the currently served model, including its title, description, identifiers, GitHub and S3 links, and Docker Hub details. A model must be served before running this command.",
    )
    @click.option(
        "--output",
        "-o",
        default=None,
        type=click.STRING,
        help="Save model information to a file. Accepted formats: .json, .csv.",
    )
    def info(output):
        from ... import ErsiliaModel
        from ...core.session import Session
        from ...hub.content.information import InformationDisplayer
        from ..messages import no_model_served, wrong_extension

        session = Session(config_json=None)
        model_id = session.current_model_id()
        service_class = session.current_service_class()
        if model_id is None:
            no_model_served()
            return
        if output and not output.endswith((".json", ".csv")):
            wrong_extension([".json", ".csv"], err=True)
            return
        mdl = ErsiliaModel(model_id, service_class=service_class)
        info = mdl.info()
        if output:
            if output.endswith(".json"):
                with open(output, "w") as f:
                    json.dump(info, f, indent=4)
            else:
                with open(output, "w", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(["Field", "Value"])
                    for key, value in info.items():
                        writer.writerow(
                            [
                                key,
                                value
                                if not isinstance(value, list)
                                else ", ".join(str(v) for v in value),
                            ]
                        )
            echo(f"Model information saved to {output}.", fg="green")
        else:
            InformationDisplayer(info).echo()

    return info
