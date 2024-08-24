from ersilia.default import INFERENCE_STORE_API_URL
import requests
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

from pydantic import BaseModel
class InferenceStoreApiPayload(BaseModel):
    model: str
    inputs: list[str] = [] # validation error if inputs is None e.g. if inchi->smiles fails


def store_has_model(model_id: str) -> bool:
    response = requests.get(
        INFERENCE_STORE_API_URL + "/model",
        params={
            "modelid": model_id
        },
        timeout=60
        )
    if response.status_code == 200:
        print("Model found in inference store: ", response.json())
        return True
    print("Model not found in inference store")
    return False
