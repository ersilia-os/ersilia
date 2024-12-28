import os
import sys

import click
import requests

from ersilia.default import INFERENCE_STORE_API_URL


class InferenceStoreMessage(object):
    """
    Base class for inference store messages.

    Parameters
    ----------
    model_id : str
        The ID of the model for which the message is being generated.
    """

    def __init__(self, model_id):
        self.model_id = model_id

    def _echo(self, text, **styles):
        return click.echo(click.style(text, **styles))


class OutputSource:
    """
    Class to define output source options.
    """

    LOCAL_ONLY = "local-only"
    CLOUD_ONLY = "cloud-only"
    ALL = [
        LOCAL_ONLY,
        CLOUD_ONLY,
    ]

    @classmethod
    def is_local(cls, option):
        """
        Check if the option is local.

        Parameters
        ----------
        option : str
            The option to check.

        Returns
        -------
        bool
            True if the option is local, False otherwise.
        """
        return option == cls.LOCAL_ONLY

    @classmethod
    def is_cloud(cls, option):
        """
        Check if the option is cloud.

        Parameters
        ----------
        option : str
            The option to check.

        Returns
        -------
        bool
            True if the option is cloud, False otherwise.
        """
        return option == cls.CLOUD_ONLY


class ModelNotInStore(InferenceStoreMessage):
    """
    Message class for models not found in the inference store.

    Parameters
    ----------
    model_id : str
        The ID of the model that is not found.
    """

    def __init__(self, model_id):
        super().__init__(model_id)
        self.model_id = model_id

    def echo(self):
        """
        Echo the message for model not found in inference store.
        """
        super()._echo(
            "Model {0} could not be found in inference store".format(self.model_id),
            fg="red",
        )
        super()._echo(
            "Please serve the model locally: ersilia serve {0} --output-source {1}".format(
                self.model_id, OutputSource.LOCAL_ONLY
            )
        )
        sys.exit(0)


class PrecalculationsNotInStore(InferenceStoreMessage):
    """
    Message class for precalculations not found in the inference store.

    Parameters
    ----------
    model_id : str
        The ID of the model for which precalculations are not found.
    """

    def __init__(self, model_id):
        super().__init__(model_id)
        self.model_id = model_id

    def echo(self):
        """
        Echo the message for precalculations not found in inference store.
        """
        super()._echo(
            "Precalculations for model {0} could not be found in inference store".format(
                self.model_id
            ),
            fg="red",
        )
        super()._echo(
            "Please serve the model locally: ersilia serve {0} --output-source {1}".format(
                self.model_id, OutputSource.LOCAL_ONLY
            )
        )
        sys.exit(0)


class PrecalculationsInStore(InferenceStoreMessage):
    """
    Message class for precalculations found in the inference store.

    Parameters
    ----------
    model_id : str
        The ID of the model for which precalculations are found.
    output_url : str
        The URL where the precalculations can be downloaded.
    """

    def __init__(self, model_id, output_url):
        super().__init__(model_id)
        self.output_url = output_url

    def echo(self):
        """
        Echo the message for precalculations available for download.

        Parameters
        ----------
        output_url : str
            The URL for downloading the precalculations.
        """
        super()._echo(
            "Precalculations for model {0} are now available for download via this link (expires in 60 minutes): {1}".format(
                self.model_id, self.output_url
            ),
            fg="green",
        )
        sys.exit(0)


def store_has_model(model_id: str) -> bool:
    """
    Check if the model exists in the inference store.

    Parameters
    ----------
    model_id : str
        The ID of the model to check.

    Returns
    -------
    bool
        True if the model exists in the store, False otherwise.
    """
    response = requests.get(
        INFERENCE_STORE_API_URL + "/model", params={"modelid": model_id}, timeout=60
    )
    if response.status_code == 200:
        print(f"Model {model_id} found in inference store")
        return True
    print(f"Model {model_id} not found in inference store")
    return False


def delete_file_upon_upload(response_code: int, file_path: str):
    """
    Delete the file upon successful upload.

    Parameters
    ----------
    response_code : int
        The HTTP response code from the upload request.
    file_path : str
        The path of the file to delete.
    """
    if response_code == 200:
        try:
            os.remove(file_path)
            print(f"File {file_path} deleted successfully.")
        except Exception as e:
            print(
                f"Failed to delete file {file_path}, please delete manually. Error: {e}"
            )
