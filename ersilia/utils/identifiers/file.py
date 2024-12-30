try:
    from Crypto.Hash import MD5
except:
    MD5 = None


class FileIdentifier(object):
    """
    A class to handle file identification and generate MD5 hashes for files.

    Parameters
    ----------
    chunk_size : int, optional
        The size of the chunks to read from the file. Default is 10000 bytes.
    """

    def __init__(self, chunk_size=10000):
        self.chunk_size = chunk_size

    def encode(self, filename, n=None):
        """
        Generate an MD5 hash for the given file.

        Parameters
        ----------
        filename : str
            The path to the file to hash.
        n : int, optional
            The number of characters of the hash to return. Default is None, which returns the full hash.

        Returns
        -------
        str
            The MD5 hash of the file, or the filename if MD5 is not available.
        """
        if MD5 is None:
            return filename
        else:
            h = MD5.new()
            with open(filename, "rb") as f:
                while True:
                    chunk = f.read(self.chunk_size)
                    if len(chunk):
                        h.update(chunk)
                    else:
                        break
            return h.hexdigest()[:n]


Identifier = FileIdentifier
