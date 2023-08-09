try:
    from Crypto.Hash import MD5
except:
    MD5 = None


class FileIdentifier(object):
    def __init__(self, chunk_size=10000):
        self.chunk_size = chunk_size

    def encode(self, filename, n=None):
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
