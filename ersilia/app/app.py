import os
import subprocess
from .. import ErsiliaBase


class AppBase(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def app_script(self, model_id):
        filename = os.path.join(self._dest_dir, model_id, "app.py")
        return filename

    def _is_swagger(self, model_id):
        filename = self.app_script(model_id)
        if not os.path.exists(filename):
            return True
        else:
            return False

    def _is_streamlit(self, model_id):
        filename = self.app_script(model_id)
        with open(filename, "r") as f:
            text = f.read()
        if "import streamlit" in text:
            return True
        else:
            return False

    def _is_dash(self, model_id):
        filename = self.app_script(model_id)
        with open(filename, "r") as f:
            text = f.read()
        if "import dash" in text:
            return True
        else:
            return False


class SwaggerApp(AppBase):

    def __init__(self, config_json=None):
        AppBase.__init__(self, config_json=config_json)

    def run(self, model_id):
        pass # TODO


class StreamlitApp(AppBase):

    def __init__(self, config_json=None):
        AppBase.__init__(self, config_json=config_json)

    def run(self, model_id):
        if not self._is_streamlit(model_id):
            return 0
        filename = os.path.join(self._dest_dir, model_id, "app.py")
        if os.path.exists(filename):
            subprocess.Popen("streamlit run %s" % filename, shell=True).wait()
            return 1
        else:
            return 0


class DashApp(AppBase):

    def __init__(self, config_json=None):
        AppBase.__init__(self, config_json=config_json)

    def run(self):
        if not self._is_dash(model_id):
            return 0
        pass # TODO
