from Crypto.Hash import MD5


class FileIdentifier(object):

    def __init__(self, chunk_size=10000):
        super().__init__()
        self.chunk_size = chunk_size

    def encode(self, filename, n=None):
        h = MD5.new()
        with open(filename, "rb") as f:
            while True:
                chunk = f.read(self.chunk_size)
                if len(chunk):
                    h.update(chunk)
                else:
                    break
        return h.hexdigest()[:n]
