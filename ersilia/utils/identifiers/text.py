import hashlib


class TextIdentifier(object):
    def __init__(self):
        pass

    def encode(self, text: str) -> str:
        return hashlib.md5(text.encode("utf-8")).hexdigest()
