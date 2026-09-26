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

    def _serving(session, mdl):
        # How the model is being served in this terminal, for the Serving
        # section: the same details as the 'Model served' summary.
        import os

        from ...default import SERVICE_CLASS_LABELS
        from ...utils import tmp_pid_file
        from ...utils.ports import normalize_connect_url

        data = session.get() or {}
        url = pid = container = None
        pid_file = tmp_pid_file(mdl.model_id)
        if os.path.isfile(pid_file):
            with open(pid_file) as f:
                lines = [line.split() for line in f if line.strip()]
            if lines:
                pid, url = lines[-1][0], lines[-1][1]
                container = lines[-1][2] if len(lines[-1]) > 2 else None
        rows = {}
        if url:
            url = normalize_connect_url(url).rstrip("/")
            rows["URL"] = f"[link={url}][cyan]{url}[/cyan][/link]"
            rows["Docs"] = f"[link={url}/docs][cyan]{url}/docs[/cyan][/link]"
        if pid and pid != "-1":
            rows["PID"] = pid
        if container and container != "-":
            rows["Container"] = container
        rows["Session"] = session._session_dir
        service = data.get("service_class")
        rows["Service"] = SERVICE_CLASS_LABELS.get(service, service)
        try:
            rows["Endpoints"] = ", ".join(mdl.get_apis())
        except Exception:
            pass
        read, write = data.get("read_store"), data.get("write_store")
        store = "Disabled"
        if read or write:
            store = "Enabled: " + " & ".join(
                n for n, on in (("Read", read), ("Write", write)) if on
            )
        rows["Store"] = store
        rows["Local cache"] = "Enabled" if data.get("local_cache") else "Disabled"
        rows["Tracking"] = "Enabled" if data.get("track_runs") else "Disabled"
        return rows

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
            no_model_served(
                hint="To read a model's card without serving it, use 'ersilia catalog --card MODEL'."
            )
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
            InformationDisplayer(info, serving=_serving(session, mdl)).echo()

    return info
