from .. import ErsiliaBase
from ..io.input import GenericInputAdapter
from .utils import InferenceStoreApiPayload
import requests, uuid

INFERENCE_STORE_UPLOAD_API_URL = ""

class InferenceStoreApi(ErsiliaBase):

    def __init__(self, model_id):
        ErsiliaBase.__init__(self)
        self.model_id = model_id
        self._upload_url = None
        self.input_adapter = GenericInputAdapter(model_id=self.model_id)

    ### Call Lambda function to get upload URL ###
    def _get_upload_url(self) -> str:
        try:
            upload_url = requests.get(
                INFERENCE_STORE_API_URL, params={"model_id": self.model_id}
            )
        except:
            upload_url = None
        return upload_url

    ### TODO: Send payload to S3 ###
    def _post_inputs(self, payload, request_id) -> str:
        return ""
    
    ### TODO: Get precalculations from S3 ###
    def _get_outputs(self, model_id, request_id) -> dict:
        return {}

    def has_model(self) -> bool:
        if self._upload_url is None: # TODO: maybe another condition here to check if URL has expired
            self._upload_url = self._get_upload_url()
        return self._upload_url is not None

    def get_precalculations(self, input: str):

        def generate_request_id():
            return str(uuid.uuid4())

        # print("-----------------------")
        # print("GOT MODEL ID")
        # print(model_id)
        # print("GOT INPUT")
        # print(input)
        # print("ADAPTED INPUT")
        adapted_input_generator = self.input_adapter.adapt_one_by_one(input)
        smiles_list = [val["input"] for val in adapted_input_generator]
        payload = InferenceStoreApiPayload(model=self.model_id, inputs=smiles_list)
        # print(payload.model_id)
        # print(payload.inputs)
        # print(payload.model_dump())
        
        request_id = generate_request_id()

        response = self._post_inputs(payload, request_id)
        if response.status_code == 200:
            print('File uploaded successfully')
        else:
            print('Failed to upload file:', response.status_code, response.text)

        precalculations = self._get_outputs(self.model_id, request_id)

        return "True"
