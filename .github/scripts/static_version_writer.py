import os
import re
from importlib.util import module_from_spec, spec_from_file_location

ROOT = os.path.dirname(os.path.abspath(__file__))


def version_writer(func):
    def wrapper(package_path):
        version = func(package_path)
        toml_path = os.path.join(package_path, "..", "pyproject.toml")
        with open(toml_path, "r") as f:
            toml_content = f.read()
        toml_content = re.sub(
            r'version\s=\s"[0-9\.]+"', f'version = "{version}"', toml_content
        )
        with open(toml_path, "w") as f:
            f.write(toml_content)
        return version

    return wrapper


@version_writer
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
