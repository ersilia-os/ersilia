"""
The Ersilia Python API.

``Model`` and ``Catalog`` mirror the CLI commands (``ersilia fetch``,
``serve``, ``run``, ``info``, ``example``, ``close``, ``delete`` and
``catalog``), with the same options, defaults and rules. As a library, they
print nothing unless ``verbose=True``, never prompt, never end the process,
and raise an ``ErsiliaError`` subclass when something goes wrong.
"""

import json
import os
import re
import shutil
import tempfile

from ..default import DEFAULT_BATCH_SIZE, INFORMATION_FILE
from ._runtime import library_call

FETCHED = "Model fetched successfully"
DOCKER_SERVICES = ("pulled_docker", "docker")
TRACKING_USE_CASES = ("local", "self-service", "hosted", "test")
EXAMPLE_MODES = ("random", "curated", "deterministic")
CATALOG_TASKS = ("Annotation", "Representation", "Sampling")
ACCESS_LEVELS = ("public", "private")


def _check_extension(path, allowed, what="output"):
    from ..utils.exceptions_utils.api_exceptions import InvalidOptionError

    if not str(path).endswith(tuple(allowed)):
        raise InvalidOptionError(
            "The {0} file must end in {1}.".format(what, " or ".join(allowed))
        )


def _check_choice(value, allowed, name):
    from ..utils.exceptions_utils.api_exceptions import InvalidOptionError

    if value is None:
        return None
    for choice in allowed:
        if str(value).lower() == choice.lower():
            return choice
    raise InvalidOptionError(
        "{0} must be one of: {1}.".format(name, ", ".join(allowed)),
        "Got {0!r}.".format(value),
    )


class Model:
    """
    An Ersilia model, identified by its identifier or slug.

    Parameters
    ----------
    model : str
        Model identifier (e.g. ``"eos4e40"``) or slug (e.g. ``"chemprop-antibiotic"``).
    verbose : bool, optional
        Show the same progress lines as the CLI. By default nothing is printed.

    Attributes
    ----------
    model_id : str
        The model identifier.
    slug : str
        The model slug.

    Raises
    ------
    InvalidModelIdentifierError
        If the model is not in the Ersilia Model Hub.

    Examples
    --------
    .. code-block:: python

        from ersilia.api import Model

        model = Model("eos4e40")
        model.fetch()
        model.serve()
        df = model.run(["CCO", "c1ccccc1"])
        model.close()

        # Or serve and close automatically:
        with Model("eos4e40") as model:
            df = model.run("input.csv")
    """

    def __init__(self, model, verbose=False):
        from ..core.modelbase import ModelBase

        self.verbose = verbose
        with library_call(verbose):
            mb = ModelBase(model)
        self.model_id = mb.model_id
        self.slug = mb.slug

    def __repr__(self):
        return "Model({0!r})".format(self.model_id)

    # Helpers

    def _session(self):
        from ..core.session import Session

        return Session(config_json=None)

    def _served_now(self):
        # The model this session serves and whether it still runs, like the
        # CLI: a record whose server stopped is cleared, as 'ersilia' does.
        session = self._session()
        served, status = session.served_model()
        if status == "stale":
            session.clear_stale(served)
            return served, False
        return served, status == "running"

    def _served_here(self):
        served, running = self._served_now()
        return running and served == self.model_id

    def _require_served(self):
        from ..utils.exceptions_utils.api_exceptions import ModelNotServedError

        served, running = self._served_now()
        if served == self.model_id and not running:
            raise ModelNotServedError(self.model_id, stopped=True)
        if not (running and served == self.model_id):
            raise ModelNotServedError(self.model_id)

    def _serving_info(self):
        # What 'ersilia info' shows in its Serving section.
        from ..utils import tmp_pid_file
        from ..utils.ports import normalize_connect_url

        session = self._session()
        data = session.get() or {}
        url = pid = container = None
        pid_file = tmp_pid_file(self.model_id)
        if os.path.isfile(pid_file):
            with open(pid_file) as f:
                lines = [line.split() for line in f if line.strip()]
            if lines:
                pid, url = lines[-1][0], lines[-1][1]
                container = lines[-1][2] if len(lines[-1]) > 2 else None
        url = normalize_connect_url(url).rstrip("/") if url else None
        try:
            apis = self._served_model().get_apis()
        except Exception:
            apis = []
        return {
            "model_id": self.model_id,
            "url": url,
            "docs": "{0}/docs".format(url) if url else None,
            "service": data.get("service_class"),
            "pid": int(pid) if pid not in (None, "-1") else None,
            "container": container if container not in (None, "-") else None,
            "session": session._session_dir,
            "apis": apis,
        }

    def _require_fetched(self):
        from ..utils.exceptions_utils.exceptions import ModelNotAvailableLocallyError

        if not self.is_fetched():
            raise ModelNotAvailableLocallyError(self.model_id)

    def _served_model(self):
        # An ErsiliaModel bound to what this session serves. It never
        # prompts or fetches: the model must be available locally.
        from ..core.model import ErsiliaModel

        session = self._session()
        return ErsiliaModel(
            self.model_id,
            output_source=session.current_output_source(),
            service_class=session.current_service_class(),
            config_json=None,
            fetch_if_not_available=False,
        )

    # Commands

    def fetch(
        self,
        from_dir=None,
        from_github=False,
        from_s3=False,
        from_hosted=None,
        version=None,
    ):
        """
        Fetch the model, like ``ersilia fetch``.

        By default the model is fetched from DockerHub. Give one of the other
        sources to fetch it from there instead.

        Parameters
        ----------
        from_dir : str, optional
            Fetch from a local copy of the model repository.
        from_github : bool, optional
            Fetch from GitHub.
        from_s3 : bool, optional
            Fetch from Ersilia's S3 bucket.
        from_hosted : str, optional
            URL of a hosted model to connect to.
        version : str, optional
            Docker image tag to fetch from DockerHub (default: the latest).

        Returns
        -------
        bool
            True if the model was fetched now, False if it was already fetched.

        Raises
        ------
        ModelFetchError
            If the model could not be fetched.
        """
        from ..hub.fetch.fetch import ALREADY_FETCHED, ModelFetcher
        from ..utils.asyncio_utils import run_coroutine
        from ..utils.echo import echo
        from ..utils.exceptions_utils.api_exceptions import (
            InvalidOptionError,
            ModelFetchError,
        )

        chosen = [
            name
            for name, on in (
                ("from_dir", from_dir),
                ("from_github", from_github),
                ("from_s3", from_s3),
                ("from_hosted", from_hosted),
            )
            if on
        ]
        if len(chosen) > 1:
            raise InvalidOptionError(
                "Choose only one source; got {0}.".format(", ".join(chosen))
            )
        from_dockerhub = not chosen
        with library_call(self.verbose):
            if version is not None and not from_dockerhub:
                echo(
                    "version only applies to DockerHub, so it is ignored.", fg="yellow"
                )
            if from_dir is not None:
                from ..utils.checks import check_fetch_folder

                check_fetch_folder(self.model_id, from_dir)
            mf = ModelFetcher(
                repo_path=from_dir,
                force_from_github=from_github,
                force_from_s3=from_s3,
                force_from_dockerhub=from_dockerhub,
                img_version=version,
                force_from_hosted=from_hosted is not None,
                hosted_url=from_hosted,
                local_dir=from_dir,
                slug=self.slug,
            )
            result = run_coroutine(mf.fetch(self.model_id))
        if result.fetch_success:
            return result.reason == FETCHED
        if result.reason == ALREADY_FETCHED:
            return False
        raise ModelFetchError(self.model_id, result.reason)

    def serve(
        self,
        port=None,
        track=False,
        tracking_use_case="local",
        enable_cache=False,
        read_store=False,
        write_store=False,
        access=None,
        nearest_neighbors=False,
        max_cache_memory_frac=None,
    ):
        """
        Serve the model, like ``ersilia serve``.

        Parameters
        ----------
        port : int, optional
            Port for the model server. By default a free port is used.
        track : bool, optional
            Track runs to monitor model and system performance.
        tracking_use_case : str, optional
            With ``track``: one of "local", "self-service", "hosted" or "test".
        enable_cache : bool, optional
            Cache results in a local Redis store.
        read_store, write_store : bool, optional
            Read results from, or write them to, the Isaura store.
        access : str, optional
            "public" or "private"; needed to write to the store.
        nearest_neighbors : bool, optional
            Use nearest-neighbor matches when reading from the store.
        max_cache_memory_frac : float, optional
            Maximum fraction (0.0-1.0) of system RAM the cache may use.

        Returns
        -------
        dict
            How the model is served, as in the Serving section of
            ``ersilia info``: ``model_id``, ``url``, ``docs``, ``service``,
            ``pid`` (None for containers), ``container``, ``session`` and
            ``apis``. If the model is already served here, it is not started
            again and these details are returned.

        Raises
        ------
        ModelNotAvailableLocallyError
            If the model is not fetched.
        SessionBusyError
            If another model is already served in this session.
        InvalidOptionError, MissingDependencyError, ModelServeError
            For invalid options, a missing Isaura installation, or a failed start.
        """
        from ..core.model import ErsiliaModel
        from ..utils.exceptions_utils.api_exceptions import (
            InvalidOptionError,
            ModelServeError,
            SessionBusyError,
        )
        from ..utils.exceptions_utils.exceptions import MissingDependencyError
        from ..utils.session import register_model_session

        tracking_use_case = _check_choice(
            tracking_use_case, TRACKING_USE_CASES, "tracking_use_case"
        )
        access = _check_choice(access, ACCESS_LEVELS, "access")
        if max_cache_memory_frac is not None and not 0 < max_cache_memory_frac <= 1:
            raise InvalidOptionError(
                "max_cache_memory_frac must be between 0.0 and 1.0.",
                "Got {0!r}.".format(max_cache_memory_frac),
            )
        if write_store and access is None:
            raise InvalidOptionError(
                "Writing to the store needs an access level.",
                "Pass access='public' or access='private'.",
            )
        if read_store or write_store:
            try:
                import isaura  # noqa: F401
            except ImportError:
                raise MissingDependencyError("isaura")
        served, running = self._served_now()
        if running and served == self.model_id:
            # Already served here: nothing to do, as 'ersilia serve' says.
            return self._serving_info()
        if running and served is not None:
            raise SessionBusyError(served, self.model_id)
        self._require_fetched()

        with library_call(self.verbose):
            mdl = ErsiliaModel(
                self.model_id,
                output_source=None,
                preferred_port=port,
                cache=enable_cache,
                maxmemory=max_cache_memory_frac,
                read_store=read_store,
                write_store=write_store,
                access=access,
                nearest_neighbors=nearest_neighbors,
                fetch_if_not_available=False,
            )
            mdl.serve(track_runs=tracking_use_case if track else None)
            if mdl.url is None:
                raise ModelServeError(self.model_id)
            register_model_session(mdl.model_id, mdl.session._session_dir)
        return self._serving_info()

    def run(self, input, output=None, batch_size=DEFAULT_BATCH_SIZE):
        """
        Run the served model, like ``ersilia run``.

        Parameters
        ----------
        input : str or list of str
            A CSV file with one column of inputs, or a list of inputs.
        output : str, optional
            A .csv or .h5 file to write the results to. If not given, the
            results are returned as a DataFrame.
        batch_size : int, optional
            Number of inputs sent to the model at a time (at least 1).

        Returns
        -------
        pandas.DataFrame or str
            The results, or the path of ``output`` when it is given.

        Raises
        ------
        ModelNotServedError
            If this model is not served in this session, or it stopped.
        InvalidOptionError
            For a problem with the input or output file (the same checks as
            ``ersilia run``: missing or empty file, several columns, wrong
            encoding or separator, missing output folder, and so on).
        EmptyRunOutputError
            If the model produced no output.
        """
        import pandas as pd

        from ..utils.exceptions_utils.api_exceptions import (
            EmptyRunOutputError,
            InvalidOptionError,
        )

        if not isinstance(batch_size, int) or batch_size < 1:
            raise InvalidOptionError(
                "batch_size must be a whole number of at least 1.",
                "Got {0!r}.".format(batch_size),
            )
        self._require_served()
        if output is not None:
            _check_extension(output, [".csv", ".h5"])
            ids = re.findall(r"eos[0-9][a-z0-9]{3}", os.path.basename(output))
            if ids and ids[0] != self.model_id:
                raise InvalidOptionError(
                    "The output file name mentions {0}, but this model is {1}.".format(
                        ids[0], self.model_id
                    )
                )
        from ..utils.checks import check_run_arguments

        tmp_dir = tempfile.mkdtemp(prefix="ersilia-api-")
        extracted = None
        try:
            if isinstance(input, (list, tuple)):
                if not input:
                    raise InvalidOptionError("No inputs were given.")
                input_path = os.path.join(tmp_dir, "input.csv")
                pd.DataFrame({"input": list(input)}).to_csv(input_path, index=False)
            elif isinstance(input, str):
                input_path = input
            else:
                raise InvalidOptionError(
                    "The input must be a CSV file with one column of inputs, or a list of inputs."
                )
            output_path = output or os.path.join(tmp_dir, "output.csv")
            with library_call(self.verbose):
                # The same input and output checks as 'ersilia run'.
                input_path, extracted = check_run_arguments(input_path, output_path)
                self._served_model().run(
                    input=input_path, output=output_path, batch_size=batch_size
                )
            if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
                raise EmptyRunOutputError(self.model_id)
            if output is not None:
                return output
            return pd.read_csv(output_path)
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)
            if extracted:
                shutil.rmtree(extracted, ignore_errors=True)

    def info(self, output=None):
        """
        Get the model information, like ``ersilia info``.

        The model must be fetched; it does not need to be served.

        Parameters
        ----------
        output : str, optional
            Also write the information to a .json or .csv file.

        Returns
        -------
        dict
            The model information.

        Raises
        ------
        ModelNotAvailableLocallyError
            If the model is not fetched.
        """
        from ..core.base import ErsiliaBase
        from ..hub.content.information import write_fields

        if output is not None:
            _check_extension(output, [".json", ".csv"])
        self._require_fetched()
        with library_call(self.verbose):
            path = os.path.join(
                ErsiliaBase()._model_path(self.model_id), INFORMATION_FILE
            )
            with open(path, "r") as f:
                info = json.load(f)
        if output is not None:
            write_fields(info, output)
        return info

    def example(self, n_samples=5, mode="random", output=None):
        """
        Generate example inputs for the model, like ``ersilia example``.

        The model must be fetched; it does not need to be served.

        Parameters
        ----------
        n_samples : int, optional
            Number of examples (ignored in "curated" mode).
        mode : str, optional
            "random", "curated" (the model's own examples) or "deterministic".
        output : str, optional
            Also write the examples to this CSV file.

        Returns
        -------
        list of str
            The example inputs.

        Raises
        ------
        ModelNotAvailableLocallyError
            If the model is not fetched.
        InvalidOptionError
            For an unknown mode, fewer than 1 example, a missing output
            folder, or curated mode on a model without curated examples.
        """
        import pandas as pd

        from ..io.input import ExampleGenerator
        from ..utils.exceptions_utils.api_exceptions import InvalidOptionError

        mode = _check_choice(mode, EXAMPLE_MODES, "mode")
        if mode != "curated" and (not isinstance(n_samples, int) or n_samples < 1):
            raise InvalidOptionError(
                "n_samples must be a whole number of at least 1.",
                "Got {0!r}.".format(n_samples),
            )
        if output is not None:
            _check_extension(output, [".csv"])
            folder = os.path.dirname(output) or "."
            if not os.path.isdir(folder):
                raise InvalidOptionError(
                    "The output folder {0} does not exist.".format(folder)
                )
        self._require_fetched()
        if mode == "curated":
            from ..utils.checks import check_curated_examples

            with library_call(self.verbose):
                check_curated_examples(self.model_id)
        tmp_dir = tempfile.mkdtemp(prefix="ersilia-api-")
        try:
            path = output or os.path.join(tmp_dir, "examples.csv")
            with library_call(self.verbose):
                ExampleGenerator(model_id=self.model_id).example(
                    n_samples, path, mode=mode
                )
            return pd.read_csv(path).iloc[:, 0].astype(str).tolist()
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    def close(self):
        """
        Close the served model, like ``ersilia close``.

        Returns
        -------
        bool
            True if the model was closed (also when it had already stopped),
            False if no model is served in this session, so there was nothing
            to close. A session named with ERSILIA_SESSION is removed.

        Raises
        ------
        ModelNotServedError
            If another model is served in this session.
        DockerNotActiveError
            If the model runs in Docker and Docker is not running.
        """
        from ..utils.exceptions_utils.api_exceptions import ModelNotServedError
        from ..utils.exceptions_utils.serve_exceptions import DockerNotActiveError
        from ..utils.session import deregister_model_session, end_named_session

        session = self._session()
        served, status = session.served_model()
        if served is None:
            # Nothing to close, as 'ersilia close' says.
            end_named_session()
            return False
        if served != self.model_id:
            raise ModelNotServedError(self.model_id)
        if status == "stale":
            # It had already stopped: forget it, as 'ersilia close' does.
            session.clear_stale(served)
            end_named_session()
            return True
        if session.current_service_class() in DOCKER_SERVICES:
            from ..setup.requirements.docker import DockerRequirement

            if not DockerRequirement().is_active():
                # Keep the record: the container stops along with Docker.
                raise DockerNotActiveError()
        with library_call(self.verbose):
            self._served_model().close()
            deregister_model_session(self.model_id)
            end_named_session()
        return True

    def delete(self):
        """
        Delete the model from this computer, like ``ersilia delete``.

        Raises
        ------
        ModelNotAvailableLocallyError
            If the model is not fetched.
        ModelNotDeletableError
            If the model cannot be deleted: it is served in this session or
            another one, or Docker is not running.
        """
        from ..hub.delete.delete import ModelFullDeleter
        from ..utils.exceptions_utils.api_exceptions import ModelNotDeletableError
        from ..utils.exceptions_utils.exceptions import ModelNotAvailableLocallyError

        if self._served_here():
            # 'ersilia delete' asks first (default no); a library never asks.
            raise ModelNotDeletableError(
                self.model_id,
                "It is being served in this session. Call close() first.",
            )
        with library_call(self.verbose):
            md = ModelFullDeleter()
            can_delete, reason = md.can_be_deleted(self.model_id)
            if not can_delete:
                if "not available locally" in reason:
                    raise ModelNotAvailableLocallyError(self.model_id)
                raise ModelNotDeletableError(self.model_id, reason)
            md.delete(self.model_id)

    def is_fetched(self):
        """
        Tell whether the model is fetched.

        Returns
        -------
        bool
            True if the model is available locally.
        """
        from ..core.modelbase import ModelBase

        with library_call(self.verbose):
            return bool(ModelBase(self.model_id).is_available_locally())

    def __enter__(self):
        if not self._served_here():
            self.serve()
        return self

    def __exit__(self, exc_type, exc, tb):
        if self._served_here():
            self.close()
        return False


class Catalog:
    """
    The Ersilia Model Hub catalog, like ``ersilia catalog``.

    Parameters
    ----------
    verbose : bool, optional
        Show the same progress lines as the CLI. By default nothing is printed.

    Examples
    --------
    .. code-block:: python

        from ersilia.api import Catalog

        catalog = Catalog()
        hub = catalog.hub(task="Annotation")
        local = catalog.local()
        card = catalog.card("eos4e40")
    """

    def __init__(self, verbose=False):
        self.verbose = verbose
        # Send Ersilia's log records to its log files from now on.
        with library_call(verbose):
            pass

    def _table(self, hub, more, task, output):
        import pandas as pd

        from ..hub.content.catalog import ModelCatalog

        task = _check_choice(task, CATALOG_TASKS, "task")
        if output is not None:
            _check_extension(output, [".csv", ".json"])
        with library_call(self.verbose):
            mc = ModelCatalog(less=not more, task=task)
            table = mc.hub() if hub else mc.local()
            if output is not None:
                table.write(output)
        return pd.DataFrame(table.data or [], columns=table.columns)

    def hub(self, more=False, task=None, output=None):
        """
        List the models in the Ersilia Model Hub, like ``ersilia catalog --hub``.

        Parameters
        ----------
        more : bool, optional
            Include more columns.
        task : str, optional
            Only "Annotation", "Representation" or "Sampling" models.
        output : str, optional
            Also write the catalog to a .csv or .json file.

        Returns
        -------
        pandas.DataFrame
            One row per model.
        """
        return self._table(True, more, task, output)

    def local(self, more=False, task=None, output=None):
        """
        List the models fetched on this computer, like ``ersilia catalog``.

        Parameters
        ----------
        more : bool, optional
            Include more columns.
        task : str, optional
            Only "Annotation", "Representation" or "Sampling" models.
        output : str, optional
            Also write the catalog to a .csv or .json file.

        Returns
        -------
        pandas.DataFrame
            One row per model; empty if no model is fetched.
        """
        return self._table(False, more, task, output)

    def card(self, model, output=None):
        """
        Get a model's card, like ``ersilia catalog --card MODEL``.

        Parameters
        ----------
        model : str
            Model identifier or slug.
        output : str, optional
            Also write the card to a .json or .csv file.

        Returns
        -------
        dict
            The model card.

        Raises
        ------
        InvalidModelIdentifierError
            If the model is not in the Ersilia Model Hub (with suggestions).
        """
        from ..hub.content.card import ModelCard
        from ..hub.content.information import write_fields
        from ..utils.exceptions_utils.exceptions import InvalidModelIdentifierError

        if output is not None:
            _check_extension(output, [".json", ".csv"])
        from ..core.modelbase import ModelBase

        with library_call(self.verbose):
            # Accepts a slug too, and suggests a close match for a typo.
            model = ModelBase(model).model_id
            card = ModelCard().get(model)
        if not card:
            raise InvalidModelIdentifierError(model)
        if output is not None:
            write_fields(card, output)
        return card
