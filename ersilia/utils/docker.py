import os
from dockerfile_parse import DockerfileParser
import tempfile
from .identifiers.long import LongIdentifier
from .terminal import run_command, run_command_check_output


class SimpleDocker(object):

    def __init__(self):
        self.identifier = LongIdentifier()

    @staticmethod
    def _image_name(org, img, tag):
        return "%s/%s:%s" % (org, img, tag)

    def images(self):
        tmp_dir = tempfile.mkdtemp()
        tmp_file = os.path.join(tmp_dir, "images.txt")
        cmd = "docker images > {0}".format(tmp_file)
        run_command(cmd, quiet=True)
        img_dict = {}
        with open(tmp_file, "r") as f:
            h = next(f)
            rep_idx = h.find("REPOSITORY")
            tag_idx = h.find("TAG")
            img_idx = h.find("IMAGE ID")
            crt_idx = h.find("CREATED")
            for l in f:
                rep = l[rep_idx:tag_idx].strip()
                tag = l[tag_idx:img_idx].strip()
                img = l[img_idx:crt_idx].strip()
                img_dict["{0}:{1}".format(rep, tag)] = img
        return img_dict

    def containers(self, only_run):
        tmp_dir = tempfile.mkdtemp()
        tmp_file = os.path.join(tmp_dir, "containers.txt")
        if not only_run:
            all_str = "-a"
        else:
            all_str = ""
        cmd = "docker ps {0} > {1}".format(all_str, tmp_file)
        run_command(cmd, quiet=True)
        cnt_dict = {}
        with open(tmp_file, "r") as f:
            h = next(f)
            cnt_idx = h.find("CONTAINER ID")
            img_idx = h.find("IMAGE")
            cmd_idx = h.find("COMMAND")
            sts_idx = h.find("STATUS")
            pts_idx = h.find("PORTS")
            nam_idx = h.find("NAMES")
            for l in f:
                cnt = l[cnt_idx:img_idx].strip()
                img = l[img_idx:cmd_idx].strip()
                sts = l[sts_idx:pts_idx].strip()
                nam = l[nam_idx:].strip()
                cnt_dict[nam] = img
        return cnt_dict

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

    def get_baseimage(self):
        return self.baseimage

    def get_runs(self):
        structure = self.structure
        runs = []
        for d in structure:
            if d["instruction"] == "RUN":
                val = d["value"]
                for v in val.split("&&"):
                    runs += [v.strip()]
        return runs
