import os
import json
import logging
import boto3
from botocore.exceptions import ClientError, NoCredentialsError
import requests

AIRTABLE_MODEL_HUB_BASE_ID = "appgxpCzCDNyGjWc8"
AIRTABLE_TABLE_ID = "tblZGe2a2XeBxrEHP"
AWS_ACCOUNT_REGION = "eu-central-1"
ERSILIA_MODEL_HUB_S3_BUCKET = "ersilia-model-hub"


def convert_airtable_to_json(
    airtable_api_key, aws_access_key_id, aws_secret_access_key
):
    headers = {"Authorization": f"Bearer {airtable_api_key}"}
    response = requests.get(
        f"https://api.airtable.com/v0/{AIRTABLE_MODEL_HUB_BASE_ID}/{AIRTABLE_TABLE_ID}",
        headers=headers,
    )

    data = response.json()
    print(f"Keys from data response: {data.keys()}")
    print(f"Offset from data response: {data.get('offset')}")
    records_models = [record["fields"] for record in data["records"]]
    models_json = json.dumps(records_models, indent=4)

    # Load JSON in AWS S3 bucket
    s3 = boto3.client(
        "s3",
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=AWS_ACCOUNT_REGION,
    )
    try:
        s3.put_object(
            Body=models_json,
            Bucket=ERSILIA_MODEL_HUB_S3_BUCKET,
            Key="models.json",
            ACL="public-read",
        )
        print("file models.json uploaded")
    except NoCredentialsError:
        logging.error("Unable to upload tracking data to AWS: Credentials not found")
    except ClientError as e:
        logging.error(e)


if __name__ == "__main__":
    print("Getting environmental variables")
    airtable_api_key = os.environ.get("AIRTABLE_API_KEY")
    aws_access_key_id = os.environ.get("AWS_ACCESS_KEY_ID")
    aws_secret_access_key = os.environ.get("AWS_SECRET_ACCESS_KEY")

    print("Converting AirTable base to JSON file")
    convert_airtable_to_json(airtable_api_key, aws_access_key_id, aws_secret_access_key)
