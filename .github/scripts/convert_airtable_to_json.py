import json
import logging
import os

import boto3
import requests
from botocore.exceptions import ClientError, NoCredentialsError

AIRTABLE_MODEL_HUB_BASE_ID = "appR6ZwgLgG8RTdoU"
AIRTABLE_TABLE_ID = "tblAfOWRbA7bI1VTB"
AWS_ACCOUNT_REGION = "eu-central-1"
ERSILIA_MODEL_HUB_S3_BUCKET = "ersilia-model-hub"


def convert_airtable_to_json(
    airtable_api_key, aws_access_key_id, aws_secret_access_key
):
    headers = {"Authorization": f"Bearer {airtable_api_key}"}
    url = (
        f"https://api.airtable.com/v0/{AIRTABLE_MODEL_HUB_BASE_ID}/{AIRTABLE_TABLE_ID}"
    )
    offset = None
    model_records = []
    
    while True:
        params={}
        if offset:
           params["offset"] = offset
          
        response = requests.get(url, headers=headers, params=params)

        if response.status_code != 200:
            print(f"Error: {response.status_code}")
            print("Response body:", response.text)
            break

        data = response.json()
        offset = data.get("offset" , None)
        model_records.extend([record["fields"] for record in data["records"]])

        if not offset:
            break

    print(f"Number of records: {len(model_records)}")
    models_json = json.dumps(model_records, indent=4)

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
        logging.error("Unable to upload data to AWS: Credentials not found")
    except ClientError as e:
        logging.error(e)


if __name__ == "__main__":
    print("Getting environmental variables")
    airtable_api_key = os.environ.get("AIRTABLE_API_KEY")
    aws_access_key_id = os.environ.get("AWS_ACCESS_KEY_ID")
    aws_secret_access_key = os.environ.get("AWS_SECRET_ACCESS_KEY")

    print("Converting AirTable base to JSON file")
    convert_airtable_to_json(airtable_api_key, aws_access_key_id, aws_secret_access_key)
