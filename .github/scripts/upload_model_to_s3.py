import sys

from ersilia.publish.s3 import S3BucketRepoUploader

model_id = sys.argv[1]
aws_access_key_id = sys.argv[2]
aws_secret_access_key = sys.argv[3]

if len(sys.argv) == 5:
    repo_path = sys.argv[4]
else:
    repo_path = None

ru = S3BucketRepoUploader(model_id=model_id)
ru.set_credentials(
    aws_access_key_id=aws_access_key_id, aws_secret_access_key=aws_secret_access_key
)
ru.upload(repo_path=repo_path)
ru.upload_zip(repo_path=repo_path)
