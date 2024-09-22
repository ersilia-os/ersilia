from ersilia.default import INFERENCE_STORE_API_URL
import requests, os

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