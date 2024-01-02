import os
import tempfile
from ..utils.docker import SimpleDocker
from ..utils.versioning import Versioner
from .. import ErsiliaBase
from .utils.clone import ErsiliaCloner


# TODO: Make sure it is used.
class SetupBaseDocker(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.docker = SimpleDocker()
        self.versions = Versioner()
        self.cloner = ErsiliaCloner(config_json=config_json)
        self.bimg = self.cfg.ENV.DOCKER.SERVER_BASE_IMAGE

    def _parse_tag(self, tag):
        if "-slim-" in tag:
            slim = True
        else:
            slim = False
        tag = tag.split("-")
        return {
            "ver": tag[0],
            "py": tag[-1],
            "python": self.versions.reformat_py(tag[-1]),
            "slim": slim,
        }

    def _get_img_name(self, org, tag):
        return "{0}/{1}:{2}".format(org, self.bimg, tag)

    def setup(self, org, tag):
        if org != "ersiliaos":
            return
        img = self.bimg
        if self.docker.exists(org, img, tag):
            return
        ptag = self._parse_tag(tag)
        # get a copy of the repository in a temporary directory
        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        tmp_repo = self.cloner.clone(tmp_folder, version=ptag["ver"])
        # write the dockerfile
        dockerfile = """
        FROM bentoml/model-server:{0}
        MAINTAINER ersilia

        ENV LC_ALL=C.UTF-8
        ENV LANG=C.UTF-8

        WORKDIR {1}

        COPY . .

        RUN pip install .
        """.format(
            tag, self.cfg.ENV.DOCKER.IMAGE_WORKDIR
        )
        path = os.path.join(tmp_repo, "Dockerfile")
        with open(path, "w") as f:
            lines = dockerfile.split("\n")
            lines = lines[1:-1]
            for l in lines:
                f.write(l[8:] + "\n")
        self.docker.build(path=tmp_repo, org=org, img=img, tag=tag)

    def delete(self, org, tag):
        pass
