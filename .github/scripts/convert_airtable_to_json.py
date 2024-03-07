import os
import json
import logging
import boto3
from botocore.exceptions import ClientError, NoCredentialsError
import requests
from ersilia.default import AIRTABLE_MODEL_HUB_BASE_ID, ERSILIA_MODEL_HUB_S3_BUCKET

AIRTABLE_TABLE_ID: 'tblZGe2a2XeBxrEHP'
AWS_ACCOUNT_REGION = "eu-central-1"

def convert_airtable_to_json(airtable_api_key, aws_access_key_id, aws_secret_access_key):
 
    headers = {'Authorization': f'Bearer {airtable_api_key}'}
    response= requests.get(f'https://api.airtable.com/v0/{AIRTABLE_MODEL_HUB_BASE_ID}/{AIRTABLE_TABLE_ID}', headers=headers)

    data=response.json()
    records_models= [record['fields'] for record in data['records']]
    models_json=json.dump(records_models)

    #Load json file in AWS S3 bucket
    s3 = boto3.client('s3',aws_access_key_id=aws_access_key_id,aws_secret_access_key=aws_secret_access_key,region_name=AWS_ACCOUNT_REGION)
    try:
        s3.response = s3.upload_file( models_json,ERSILIA_MODEL_HUB_S3_BUCKET,'models.json',ExtraArgs={'ACL': 'public-read'})
        print("file models.json uploaded")
    except NoCredentialsError:
            logging.error("Unable to upload tracking data to AWS: Credentials not found")
    except ClientError as e:
        logging.error(e)
        return False
    
if __name__ == "__main__":

    print("Getting environmental variables")
    airtable_api_key = os.environ.get('AIRTABLE_API_KEY')
    aws_access_key_id = os.environ.get('AWS_ACCESS_KEY_ID')
    aws_secret_access_key = os.environ.get('AWS_SECRET_ACCESS_KEY')

    print("Converting AirTable base to JSON file")
    convert_airtable_to_json(airtable_api_key, aws_access_key_id, aws_secret_access_key)
