from ersilia.core.base import ErsiliaBase
from ersilia.io.input import GenericInputAdapter
from ersilia.store.utils import InferenceStoreApiPayload, store_has_model
from ersilia.default import INFERENCE_STORE_API_URL
import requests
import tempfile
import uuid


class InferenceStoreApi(ErsiliaBase):

    def __init__(self, model_id):
        ErsiliaBase.__init__(self)
        self.model_id = model_id
        self.request_id = None
        self.input_adapter = GenericInputAdapter(model_id=self.model_id)
        
    def _generate_request_id(self) -> str:
        self.request_id = str(uuid.uuid4())
        return self.request_id

    def _get_presigned_url(self) -> str:
        response = requests.get(
                INFERENCE_STORE_API_URL + "/upload-destination",
                params={
                    "modelid": self.model_id,
                    "requestid": self.request_id
                },
                timeout=60
                )
        presigned_url_response = response.json()
        return presigned_url_response

    def _post_inputs(self, inputs, presigned_url_response) -> str:
        adapted_input_generator = self.input_adapter.adapt_one_by_one(inputs)
        smiles_list = [val["input"] for val in adapted_input_generator]
        payload = InferenceStoreApiPayload(model=self.model_id, inputs=smiles_list)

        with tempfile.NamedTemporaryFile() as input_file:
            input_file.write(payload.model_dump())
            inputs = input_file.name

        with open(inputs, 'rb') as f:
            files = {'file': (inputs, f)}
            presigned_url = presigned_url_response.get('url')
            presigned_url_data = presigned_url_response.get('fields')
            http_response = requests.post(
                presigned_url,
                data=presigned_url_data,
                files=files,
                timeout=60
                )

        return http_response

    def _get_outputs(self) -> str:
        response = requests.post(
            INFERENCE_STORE_API_URL + "/precalculations",
            params={
                "requestid": self.request_id
            },
            timeout=60
            )
        return response.text

    def get_precalculations(self, inputs):
        if store_has_model(self.model_id):
            self._generate_request_id()
            presigned_url = self._get_presigned_url()

        response = self._post_inputs(inputs, presigned_url)
        if response.status_code == 204:
            print('File uploaded successfully')
        else:
            return f'Failed to upload file: {response.status_code} error ({response.text})'

        precalculations = self._get_outputs()

        # TODO: should missing keys be returned too in a file/message?
        return precalculations # this is the result returned to the CLI
