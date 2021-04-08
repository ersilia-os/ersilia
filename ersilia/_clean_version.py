import os
from setuptools.command.build_py import build_py as build_py_orig
from setuptools.command.sdist import sdist as sdist_orig

STATIC_VERSION_FILE = "_clean_static_version.py"

# No public API
__all__ = []

package_root = os.path.dirname(os.path.realpath(__file__))
package_name = os.path.basename(package_root)
distr_root = os.path.dirname(package_root)
# If the package is inside a "src" directory the
# distribution root is 1 level up.
if os.path.split(distr_root)[1] == "src":
    _package_root_inside_src = True
    distr_root = os.path.dirname(distr_root)
else:
    _package_root_inside_src = False

script_path = os.path.dirname(os.path.abspath(__file__))
clean_version_file = os.path.join(script_path, STATIC_VERSION_FILE)
with open(clean_version_file, "r") as f:
    version = f.read().split('"')[1]

__version__ = version


def _write_version(fname):
    # This could be a hard link, so try to delete it first.  Is there any way
    # to do this atomically together with opening?
    try:
        os.remove(fname)
    except OSError:
        pass
    with open(fname, "w") as f:
        f.write(
            "# This file has been created by setup.py.\n"
            "version = '{}'\n".format(__version__)
        )


class _build_py(build_py_orig):
    def run(self):
        super().run()
        _write_version(os.path.join(self.build_lib, package_name, STATIC_VERSION_FILE))


class _sdist(sdist_orig):
    def make_release_tree(self, base_dir, files):
        super().make_release_tree(base_dir, files)
        if _package_root_inside_src:
            p = os.path.join("src", package_name)
        else:
            p = package_name
        _write_version(os.path.join(base_dir, p, STATIC_VERSION_FILE))


cmdclass = dict(sdist=_sdist, build_py=_build_py)
