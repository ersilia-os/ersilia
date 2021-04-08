import zipfile
import os
import shutil


class Zipper(object):
    def __init__(self, remove):
        self.remove = remove

    def unzip(self, file, destination):
        with zipfile.ZipFile(file, "r") as zip_ref:
            zip_ref.extractall(destination)
        if self.remove:
            os.remove(file)

    def zip(self, dir_name, file):
        shutil.make_archive(file, "zip", dir_name)
        if self.remove:
            shutil.rmtree(dir_name)
