import os
import sys

import click
import requests

from ersilia.default import INFERENCE_STORE_API_URL


class InferenceStoreMessage(object):
    def __init__(self, model_id):
        self.model_id = model_id

    def _echo(self, text, **styles):
        return click.echo(click.style(text, **styles))

class OutputSource():
    LOCAL_ONLY = "local-only"
    CLOUD_ONLY = "cloud-only"
    ALL = [
        LOCAL_ONLY,
        CLOUD_ONLY,
        ]

    @classmethod
    def is_local(cls, option):
        return option == cls.LOCAL_ONLY

    @classmethod
    def is_cloud(cls, option):
        return option == cls.CLOUD_ONLY

class ModelNotInStore(InferenceStoreMessage):
    def __init__(self, model_id):
        super().__init__(model_id)
        self.model_id = model_id

    def echo(self):
        super()._echo(
            "Model {0} could not be found in inference store".format(self.model_id),
            fg="red",
        )
        super()._echo(
            "Please serve the model locally: ersilia serve {0} --output-source {1}".format(
                self.model_id,
                OutputSource.LOCAL_ONLY
                )
        )
        sys.exit(0)


class PrecalculationsNotInStore(InferenceStoreMessage):
    def __init__(self, model_id):
        super().__init__(model_id)
        self.model_id = model_id

    def echo(self):
        super()._echo(
            "Precalculations for model {0} could not be found in inference store".format(self.model_id),
            fg="red",
        )
        super()._echo(
            "Please serve the model locally: ersilia serve {0} --output-source {1}".format(
                self.model_id,
                OutputSource.LOCAL_ONLY
                )
        )
        sys.exit(0)


class PrecalculationsInStore(InferenceStoreMessage):
    def __init__(self, model_id, output_url):
        super().__init__(model_id)
        self.output_url = output_url

    def echo(self):
        super()._echo(
            "Precalculations for model {0} are now available for download via this link (expires in 60 minutes): {1}".format(
                self.model_id,
                self.output_url
                ),
            fg="green"
        )
        sys.exit(0)


def store_has_model(model_id: str) -> bool:
    response = requests.get(
        INFERENCE_STORE_API_URL + "/model",
        params={
            "modelid": model_id
        },
        timeout=60
        )
    if response.status_code == 200:
        print(f"Model {model_id} found in inference store")
        return True
    print(f"Model {model_id} not found in inference store")
    return False


def delete_file_upon_upload(response_code: int, file_path: str):
    if response_code == 200:
        try:
            os.remove(file_path)
            print(f"File {file_path} deleted successfully.")
        except Exception as e:
            print(f"Failed to delete file {file_path}, please delete manually. Error: {e}")