import os
import os.path
import shutil
from typing import Tuple

from ... import ErsiliaBase
from ...db.disk.fetched import FetchedModelsManager
from ...db.environments.localdb import EnvironmentDb
from ...db.environments.managers import DockerManager
from ...db.hubdata.localslugs import SlugDb
from ...setup.requirements.bentoml_requirement import BentoMLRequirement
from ...utils.conda import SimpleConda
from ...utils.docker import SimpleDocker
from ...utils.echo import echo  # <--- you already have semantic echo()
from ...utils.environment import Environment
from ...utils.session import (
    deregister_model_session,
    get_model_session,
    remove_session_dir,
)
from ...utils.system import is_inside_docker
from ...utils.terminal import run_command
from ..bundle.status import ModelStatus
from ..content.card import ModelCard
from ..content.catalog import ModelCatalog


def rmtree(path):
    """
    Recursively delete a directory tree.

    Parameters
    ----------
    path : str
        Path to the directory to be deleted.
    """
    try:
        shutil.rmtree(path)
    except:
        pass
    try:
        os.unlink(path)
    except:
        pass


class ModelEosDeleter(ErsiliaBase):
    """
    Deletes model data from the EOS storage directory, a directory to store fetched models locally.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the deleter.

    Methods
    -------
    delete(model_id)
        Deletes the model data from the EOS storage directory, a directory to store fetched models locally.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _model_path(self, model_id):
        folder = os.path.join(self._dest_dir, model_id)
        return folder

    def delete(self, model_id: str):
        """
        Deletes the model data from the EOS storage directory, a directory to store fetched models locally..

        Parameters
        ----------
        model_id : str
            Identifier of the model to be deleted.
        """
        folder = self._model_path(model_id)
        if not os.path.exists(folder):
            return
        self.logger.info("Removing EOS folder {0}".format(folder))
        echo(f"Removing EOS folder {folder}", fg="yellow")
        rmtree(folder)


class ModelTmpDeleter(ErsiliaBase):
    """
    Deletes temporary model data.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the deleter.

    Methods
    -------
    delete(model_id)
        Deletes the temporary model data.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _model_path(self, model_id):
        folder = os.path.join(self._tmp_dir, model_id)
        return folder

    def delete(self, model_id: str):
        """
        Deletes the temporary model data from the EOS directory.

        Parameters
        ----------
        model_id : str
            Identifier of the model to be deleted.
        """
        self.logger.debug("Attempting temporary folder delete")
        folder = self._model_path(model_id)
        if not os.path.exists(folder):
            return
        self.logger.info("Removing temporary folder {0}".format(folder))
        echo(f"Removing temporary folder {folder}", fg="yellow")
        shutil.rmtree(folder)


class ModelBundleDeleter(ErsiliaBase):
    """
    Deletes model bundles from the EOS directory.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the deleter.

    Methods
    -------
    delete(model_id)
        Deletes the model bundle.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _model_path(self, model_id):
        folder = os.path.join(self._bundles_dir, model_id)
        return folder

    def delete(self, model_id: str):
        """
        Deletes the model bundle.

        Parameters
        ----------
        model_id : str
            Identifier of the model to be deleted.
        """
        folder = self._model_path(model_id)
        if not os.path.exists(folder):
            return
        bento_folder = self._get_bentoml_location(model_id)
        if bento_folder is not None:
            self.logger.info("Removing bento folder first {0}".format(bento_folder))
            echo(f"Removing bento folder {bento_folder}", fg="yellow")
            rmtree(bento_folder)
            os.makedirs(bento_folder, exist_ok=True)
        self.logger.info("Removing bundle folder {0}".format(folder))
        echo(f"Removing bundle folder {folder}", fg="yellow")
        rmtree(folder)
        self.logger.debug("Folder removed")


class ModelBentoDeleter(ErsiliaBase):
    """
    Deletes BentoML services related to a model.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the deleter.

    Methods
    -------
    delete(model_id)
        Deletes all BentoML services related to the model.
    clean(model_id)
        Deletes all but the latest BentoML service related to the model.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _delete_service(self, service):
        cmd = "echo yes | bentoml delete %s" % service
        self.logger.debug(cmd)
        run_command(cmd)

    def _delete(self, model_id, keep_latest=True):
        bentomlreq = BentoMLRequirement()
        if bentomlreq.is_installed():
            ml = ModelCatalog()
            try:
                self.logger.debug("Looking for BentoML catalog")
                catalog = ml.bentoml()
                self.logger.debug("Catalog found")
            except:
                self.logger.debug("No BentoML Catalog available")
                catalog = None
            if catalog is None:
                return
            if len(catalog.data) == 0:
                return
            services = [r[1] for r in catalog.data if r[0] == model_id]
            if keep_latest and len(services) > 1:
                services = services[1:]
            for service in services:
                self.logger.info("Removing BentoML service {0}".format(service))
                echo(f"Removing BentoML service {service}", fg="yellow")
                self._delete_service(service)
        else:
            self.logger.debug(
                "BentoML is not installed, no need for deleting using BentoML"
            )

    def delete(self, model_id: str):
        """
        Deletes all BentoML services related to the model.

        Parameters
        ----------
        model_id : str
            Identifier of the model to be deleted.
        """
        self.logger.debug("Attempting Bento delete")
        self._delete(model_id, keep_latest=False)

    def clean(self, model_id: str):
        """
        Deletes all but the latest BentoML service related to the model.

        Parameters
        ----------
        model_id : str
            Identifier of the model to be cleaned.
        """
        self._delete(model_id, keep_latest=True)


class ModelSlugDeleter(ErsiliaBase):
    """
    Deletes model slugs from the database.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the deleter.

    Methods
    -------
    delete(model_id)
        Deletes the model slug from the database.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.slugdb = SlugDb(config_json=config_json)

    def delete(self, model_id: str):
        """
        Deletes the model slug from the database.

        Parameters
        ----------
        model_id : str
            Identifier of the model to be deleted.
        """
        echo(f"Removing slug DB entry for {model_id}", fg="yellow")
        self.slugdb.delete_by_model_id(model_id)


class ModelCondaDeleter(ErsiliaBase):
    """
    Deletes Conda environments related to a model.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the deleter.

    Methods
    -------
    delete(model_id)
        Deletes the Conda environment related to the model.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.envdb = EnvironmentDb(config_json=config_json)
        self.envdb.table = "conda"

    def _to_delete(self, model_id):
        env = self.envdb.envs_of_model(model_id)
        if len(env) != 1:
            return None
        else:
            return list(env)[0]

    def delete(self, model_id: str):
        """
        Deletes the Conda environment related to the model.

        Parameters
        ----------
        model_id : str
            Identifier of the model to be deleted.
        """
        envs = self.envdb.envs_of_model(model_id)
        envs = list(set(list(envs) + [model_id]))
        for env in envs:
            models = self.envdb.models_of_env(env)
            models.update([model_id])
            if len(models) == 1:
                self.logger.info("Deleting conda environment {0}".format(env))
                echo(f"Deleting conda environment {env}", fg="yellow")
                try:
                    conda = SimpleConda()
                    conda.delete(env)
                except:
                    continue
            self.envdb.delete(model_id, env)


class ModelPipDeleter(ErsiliaBase):
    """
    Uninstalls pip packages related to a model.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the deleter.

    Methods
    -------
    delete(model_id)
        Uninstalls the pip package related to the model.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)

    def pip_uninstall(self, model_id):
        """
        Uninstalls a Python package using pip.
        Parameters
        ----------
        model_id : str
            The name of the package to uninstall.
        """
        run_command("echo y | pip uninstall %s" % model_id)

    def delete(self, model_id: str):
        """
        Uninstalls the pip package related to the model.

        Parameters
        ----------
        model_id : str
            Identifier of the model to be deleted.
        """
        env = Environment()
        if env.has_module(model_id):
            self.logger.info("Uninstalling pip package {0}".format(model_id))
            echo(f"Uninstalling pip package {model_id}", fg="yellow")
            self.pip_uninstall(model_id)


class ModelDockerDeleter(ErsiliaBase):
    """
    Deletes Docker images and stops containers related to a model.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the deleter.

    Methods
    -------
    delete(model_id)
        Deletes Docker images and stops containers related to the model.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)

    def delete(self, model_id: str):
        """
        Deletes Docker images and stops containers related to the model.

        Parameters
        ----------
        model_id : str
            Identifier of the model to be deleted.
        """
        if is_inside_docker():
            return
        self.logger.info(
            "Removing docker images and stopping containers related to {0}".format(
                model_id
            )
        )
        echo(f"Removing docker images and containers for {model_id}", fg="yellow")
        dm = DockerManager(config_json=self.config_json)
        if dm.is_active():
            dm.delete_images(model_id)


class ModelFetchedEntryDeleter(ErsiliaBase):
    """
    Deletes fetched model entries from the database.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the deleter.

    Methods
    -------
    delete(model_id)
        Deletes the fetched model entry from the database.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.fmm = FetchedModelsManager(config_json=config_json)

    def delete(self, model_id: str):
        """
        Deletes the fetched model entry from the database.

        Parameters
        ----------
        model_id : str
            Identifier of the model to be deleted.
        """
        echo(f"Removing fetched DB entry for {model_id}", fg="yellow")
        self.fmm.delete(model_id)


class TmpCleaner(ErsiliaBase):
    """
    Cleans temporary directories.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the cleaner.

    Methods
    -------
    delete()
        Deletes the temporary directory.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def delete(self):
        """
        Deletes the temporary directory.
        """
        if os.path.exists(self._tmp_dir):
            echo(f"Cleaning tmp directory {self._tmp_dir}", fg="yellow")
            os.rmdir(self._tmp_dir)
            os.makedirs(self._tmp_dir)


class BruteDeleter(ErsiliaBase):
    """
    Brute deletion of model, once all individual deleters have been run.

    Parameters
    ----------
    model_id : str
        Identifier of the model to be deleted.
    config_json : dict, optional
        Configuration settings.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.logger.debug("Attempting brute delete")
        self.conda = SimpleConda(config_json=config_json)
        self.docker = SimpleDocker()

    def delete(self, model_id):
        """
        Deletes all remaining data related to the model.
        """
        self.model_id = model_id
        dest_dir = self._model_path(self.model_id)
        if os.path.exists(dest_dir):
            self.logger.info("Removing dest folder {0}".format(dest_dir))
            echo(f"Removing dest folder {dest_dir}", fg="yellow")
            rmtree(dest_dir)
        bundles_dir = os.path.join(self._bundles_dir, self.model_id)
        if os.path.exists(bundles_dir):
            self.logger.info("Removing bundles folder {0}".format(bundles_dir))
            echo(f"Removing bundles folder {bundles_dir}", fg="yellow")
            rmtree(bundles_dir)
        tmp_dir = os.path.join(self._tmp_dir, self.model_id)
        if os.path.exists(tmp_dir):
            self.logger.info("Removing tmp folder {0}".format(tmp_dir))
            echo(f"Removing tmp folder {tmp_dir}", fg="yellow")
            rmtree(tmp_dir)
        conda_model_ids = self.conda.list_eos_environments()
        if self.model_id in conda_model_ids:
            self.logger.info("Removing conda environment {0}".format(self.model_id))
            echo(f"Removing conda environment {self.model_id}", fg="yellow")
            self.conda.delete(self.model_id)
        docker_model_ids = self.docker.list_eos_images()
        if self.model_id in docker_model_ids:
            self.logger.info(
                "Removing docker containers related to {0}".format(self.model_id)
            )
            echo(f"Removing docker containers for {self.model_id}", fg="yellow")
            self.docker.brute_delete_by_model_id(self.model_id)


class ModelFullDeleter(ErsiliaBase):
    """
    Deletes all data related to a model.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the deleter.
    overwrite : bool, optional
        Whether to overwrite existing data, by default True.

    Methods
    -------
    needs_delete(model_id)
        Checks if the model needs to be deleted.
    delete(model_id)
        Deletes all data related to the model.
    """

    def __init__(self, config_json=None, overwrite=True):
        self.overwrite = overwrite
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.conda = SimpleConda(config_json=config_json)
        self.docker = SimpleDocker()

    def _needs_delete(self, model_id: str) -> bool:
        """
        Checks if the model needs to be deleted.

        Parameters
        ----------
        model_id : str
            Identifier of the model.

        Returns
        -------
        bool
            True if the model needs to be deleted, False otherwise.
        """
        ms = ModelStatus().status(model_id)
        for k, v in ms.items():
            if v:
                return True
        if os.path.exists(self._model_path(model_id)):
            return True
        return False

    def can_be_deleted(self, model_id: str) -> Tuple[bool, str]:
        """
        Checks if the model can be deleted.

        Parameters
        ----------
        model_id : str
            Identifier of the model.

        Returns
        -------
        bool
            True if the model can be deleted, False otherwise.
        """
        mdl_session = get_model_session(model_id)
        if mdl_session:
            self.logger.debug("Removing session {0}".format(mdl_session))
            remove_session_dir(mdl_session)
            deregister_model_session(model_id)
        needs_delete = self._needs_delete(model_id)
        mc = ModelCard(config_json=self.config_json).get(model_id)
        model_source = ModelCatalog(config_json=self.config_json)._get_model_source(mc)
        dm = DockerManager(config_json=self.config_json)
        model_ids_in_conda = self.conda.list_eos_environments()
        if model_id in model_ids_in_conda:
            return (
                True,
                "Model has a conda environment, can be deleted.",
            )
        model_ids_in_docker = self.docker.list_eos_images()
        if model_id in model_ids_in_docker:
            return (
                True,
                "Model has a docker image, can be deleted.",
            )
        if needs_delete:
            if model_source == "DockerHub" and not dm.is_active():
                return (
                    False,
                    "Model fetched through Docker but Docker engine is inactive.",
                )
            return True, "Model can be deleted."
        else:
            return (
                False,
                f"Model {model_id} is not available locally, no delete necessary.",
            )

    def delete(self, model_id: str):
        """
        Deletes all data related to the model.

        Parameters
        ----------
        model_id : str
            Identifier of the model to be deleted.
        """
        self.logger.info("Starting delete of model {0}".format(model_id))
        echo(f"Deleting all data for model {model_id}", fg="cyan")

        ModelEosDeleter(self.config_json).delete(model_id)
        ModelSlugDeleter(self.config_json).delete(model_id)
        ModelBundleDeleter(self.config_json).delete(model_id)
        ModelBentoDeleter(self.config_json).delete(model_id)
        if self.overwrite:
            ModelCondaDeleter(self.config_json).delete(model_id)
        ModelTmpDeleter(self.config_json).delete(model_id)
        ModelPipDeleter(self.config_json).delete(model_id)
        ModelDockerDeleter(self.config_json).delete(model_id)
        ModelFetchedEntryDeleter(self.config_json).delete(model_id)
        BruteDeleter(self.config_json).delete(model_id)

        self.logger.success("Model {0} deleted successfully".format(model_id))
        echo(f"Model {model_id} deleted successfully", fg="green")
