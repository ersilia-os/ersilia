import os
from dockerfile_parse import DockerfileParser
import tempfile
from .identifiers import LongIdentifier
from .terminal import run_command, run_command_check_output


class SimpleDocker(object):

    def __init__(self):
        self.identifier = LongIdentifier()

    @staticmethod
    def _image_name(org, img, tag):
        return "%s/%s:%s" % (org, img, tag)

    def exists(self, org, img, tag):
        bash_script = """
        #!/bin/bash
        image_and_tag="$1"
        image_and_tag_array=(${image_and_tag//:/ })
        if [[ "$(docker images ${image_and_tag_array[0]} | grep ${image_and_tag_array[1]} 2> /dev/null)" != "" ]]; then
          echo "True"
        else
          echo "False"
        fi
        """
        tmp_folder = tempfile.mkdtemp()
        tmp_script = os.path.join(tmp_folder, "exists.sh")
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        cmd = ["bash", tmp_script, self._image_name(org, img, tag)]
        res = run_command_check_output(cmd)
        res = res.strip()
        if res == b"False":
            return False
        if res == b"True":
            return True
        return None

    def build(self, path, org, img, tag):
        path = os.path.abspath(path)
        cwd = os.getcwd()
        os.chdir(path)
        cmd = "docker build -t %s %s" % (self._image_name(org, img, tag), path)
        run_command(cmd, quiet=True)
        os.chdir(cwd)

    def delete(self, org, img, tag):
        cmd = "docker rmi %s" % self._image_name(org, img, tag)
        run_command(cmd, quiet=True)

    def run(self, org, img, tag, name):
        if name is None:
            name = self.identifier.encode()
        cmd = "docker run -it -d --name %s %s bash" % (name, self._image_name(org, img, tag))
        run_command(cmd, quiet=True)
        return name

    @staticmethod
    def kill(name):
        cmd = "docker kill %s" % name
        run_command(cmd, quiet=True)

    @staticmethod
    def cp_from_container(name, img_path, local_path):
        local_path = os.path.abspath(local_path)
        cmd = "docker cp %s:%s %s" % (name, img_path, local_path)
        run_command(cmd, quiet=True)

    def cp_from_image(self, img_path, local_path, org, img, tag):
        name = self.run(org, img, tag, name=None)
        self.cp_from_container(name, img_path, local_path)
        self.kill(name)

    @staticmethod
    def exec_container(name, cmd):
        cmd = 'docker exec -i %s bash -c "%s"' % (name, cmd)
        run_command(cmd, quiet=True)

    def exec(self, cmd, org, img, tag, name):
        name = self.run(org, img, tag, name=name)
        self.exec_container(name, cmd)
        self.kill(name)


class SimpleDockerfileParser(DockerfileParser):

    def __init__(self, path):
        if os.path.isdir(path):
            path = os.path.join(path, "Dockerfile")
        DockerfileParser.__init__(self, path=path)

    def get_runs(self):
        structure = self.structure
        runs = []
        for d in structure:
            if d["instruction"] == "RUN":
                val = d["value"]
                for v in val.split("&&"):
                    runs += [v.strip()]
        return runs
