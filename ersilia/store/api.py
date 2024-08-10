from ersilia.core.base import ErsiliaBase
from ersilia.io.input import GenericInputAdapter
from ersilia.store.utils import InferenceStoreApiPayload
import requests, uuid

INFERENCE_STORE_UPLOAD_API_URL = ""

class InferenceStoreApi(ErsiliaBase):

    def __init__(self, model_id):
        ErsiliaBase.__init__(self)
        self.model_id = model_id
        self._upload_url = None
        self.input_adapter = GenericInputAdapter(model_id=self.model_id)

    ### Call Lambda function to get upload URL ###
    def _get_presigned_url(self) -> str:
        presigned_url = requests.get(
                INFERENCE_STORE_API_URL, params={"model_id": self.model_id}
            )
        
        # try:
        #     upload_url = requests.get(
        #         INFERENCE_STORE_API_URL, params={"model_id": self.model_id}
        #     )
        # except:
        #     upload_url = None
        # return upload_url

    def has_model(self) -> bool:
        # GET request to check if model exists

        # if self._upload_url is None: # TODO: maybe another condition here to check if URL has expired
        #     self._upload_url = self._get_upload_url()
        # return self._upload_url is not Non

    ### TODO: Send payload to S3 ###
    def _post_inputs(self, presigned_url) -> str:
        adapted_input_generator = self.input_adapter.adapt_one_by_one(input)
        smiles_list = [val["input"] for val in adapted_input_generator]
        payload = InferenceStoreApiPayload(model=self.model_id, inputs=smiles_list)
        return ""

    ### TODO: Get precalculations from S3 ###
    def _get_outputs(self, model_id=None, request_id=None) -> dict:
        return {}

    def get_precalculations(self, input: str):

        def generate_request_id():
            return str(uuid.uuid4())

        # print("-----------------------")
        # print("GOT MODEL ID")
        # print(model_id)
        # print("GOT INPUT")
        # print(input)
        # print("ADAPTED INPUT")
        
        # print(payload.model_id)
        # print(payload.inputs)
        # print(payload.model_dump())
        
        request_id = generate_request_id()
        
        presigned_url = self._get_presigned_url()

        response = self._post_inputs(presigned_url)
        if response.status_code == 200:
            print('File uploaded successfully')
        else:
            print('Failed to upload file:', response.status_code, response.text)

        precalculations = self._get_outputs(self.model_id, request_id)

        return precalculations

if __name__ == "__main__":
    abc = 'https://bzr6zxw1k0.execute-api.ap-southeast-2.amazonaws.com/4-8-24'
    response = requests.get(abc, params={"model_id": 'ey1829ey', "request_id": '12e1-2e-12e12e'})
    print(response.json())