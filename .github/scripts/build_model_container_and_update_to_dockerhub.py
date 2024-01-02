import sys
from ersilia.publish.dockerhub import DockerHubUploader


model_id = sys.argv[1]
docker_user = sys.argv[2]
docker_pwd = sys.argv[3]
du = DockerHubUploader(model_id=model_id)
du.set_credentials(docker_user=docker_user, docker_pwd=docker_pwd)
du.upload()
