import os
import os
from importlib.util import module_from_spec, spec_from_file_location

ROOT = os.path.dirname(os.path.abspath(__file__))


def get_version(package_path):
    spec = spec_from_file_location("version", os.path.join(package_path, "_version.py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    version = module.get_version_for_setup()
    return version


try:
    version = get_version(os.path.join(ROOT, "..", "..", "ersilia"))
except:
    version = get_version(os.path.join(ROOT, "ersilia"))

print(version)
