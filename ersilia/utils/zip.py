import os
import shutil
import zipfile


class Zipper(object):
    """
    A class to handle zipping and unzipping files.

    Parameters
    ----------
    remove : bool
        Whether to remove the original files after zipping/unzipping.

    Methods
    -------
    unzip(file, destination)
        Unzip a file to the specified destination.
    zip(dir_name, file)
        Zip a directory into a file.
    """

    def __init__(self, remove):
        self.remove = remove

    def unzip(self, file, destination):
        """
        Unzip a file to the specified destination.

        Parameters
        ----------
        file : str
            The path to the zip file.
        destination : str
            The destination directory.

        Returns
        -------
        None
        """
        with zipfile.ZipFile(file, "r") as zip_ref:
            zip_ref.extractall(destination)
        if self.remove:
            os.remove(file)

    def zip(self, dir_name, file):
        """
        Zip a directory into a file.

        Parameters
        ----------
        dir_name : str
            The path to the directory to zip.
        file : str
            The path to the output zip file.

        Returns
        -------
        None
        """
        shutil.make_archive(file, "zip", dir_name)
        if self.remove:
            shutil.rmtree(dir_name)
