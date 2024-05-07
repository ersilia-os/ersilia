import json
import boto3

from ...default import ERSILIA_MODEL_HUB_S3_BUCKET

AWS_ACCOUNT_REGION = "eu-central-1"


class JsonModelsInterface:
    def __init__(self, json_file_name):

        self.json_file_name = json_file_name
        self.s3_client = boto3.client("s3", region_name=AWS_ACCOUNT_REGION)

    def _read_json_file(self):

        s3_response = self.s3_client.get_object(
            Bucket=ERSILIA_MODEL_HUB_S3_BUCKET, Key=self.json_file_name
        )
        s3_object_body = s3_response.get("Body")
        content = s3_object_body.read().decode("utf-8")
        json_dict = json.loads(content)

        return json_dict

    def items(self):
        json_dict = self._read_json_file()
        for key, value in json_dict.items:
            yield key, value

    def items_all(self):
        items = self._read_json_file()
        return items
