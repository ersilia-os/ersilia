from ersilia.core.base import ErsiliaBase
from ersilia.io.input import GenericInputAdapter
from ersilia.store.utils import PrecalculationsNotInStore, PrecalculationsInStore, delete_file_upon_upload
from ersilia.default import INFERENCE_STORE_API_URL
import requests
import uuid
import csv


class InferenceStoreApi(ErsiliaBase):

    def __init__(self, model_id):
        ErsiliaBase.__init__(self)
        self.model_id = model_id
        self.request_id = None
        self.input_adapter = GenericInputAdapter(model_id=self.model_id)
        
    def _generate_request_id(self) -> str:
        self.request_id = str(uuid.uuid4())

    def _get_presigned_url(self) -> str:
        response = requests.get(
            INFERENCE_STORE_API_URL + "/upload-destination",
            params={
                "modelid": self.model_id,
                "requestid": self.request_id
            }
        )
        presigned_url_response = response.json()
        return presigned_url_response

    def _post_inputs(self, inputs, presigned_url_response) -> str:
        adapted_input_generator = self.input_adapter.adapt_one_by_one(inputs)
        smiles_list = [val["input"] for val in adapted_input_generator]

        # TODO: avoid creating a file
        ersilia_input_file = "ersilia_cloud_inputs.csv"
        with open(ersilia_input_file, "w", newline="") as input_file:
            writer = csv.writer(input_file)
            for smiles in smiles_list:
                writer.writerow([smiles])
        
        with open(ersilia_input_file, "rb") as input_file:
            files = {'file': (ersilia_input_file, input_file)}
            presigned_url = presigned_url_response.get('url')
            presigned_url_data = presigned_url_response.get('fields')
            http_response = requests.post(
                presigned_url,
                data=presigned_url_data,
                files=files
                )
        delete_file_upon_upload(http_response.status_code, ersilia_input_file)
        return http_response

    def _get_outputs(self) -> str:
        response = requests.post(
            INFERENCE_STORE_API_URL + "/predictions",
            params={
                "requestid": self.request_id,
                "modelid": self.model_id
            },
            timeout=60
            )
        output_presigned_url = response.text
        return output_presigned_url

    def get_precalculations(self, inputs):
        self._generate_request_id()
        presigned_url = self._get_presigned_url()

        response = self._post_inputs(inputs, presigned_url)
        if response.status_code == 204:
            print('File uploaded successfully')
        else:
            return f'Failed to upload file: {response.status_code} error ({response.text})'

        output_presigned_url = self._get_outputs()
        
        if not output_presigned_url:
            PrecalculationsNotInStore(self.model_id).echo()

        # this is the result returned to the CLI (should missing keys be returned too in a file/message)?
        return PrecalculationsInStore(self.model_id, output_presigned_url).echo()
