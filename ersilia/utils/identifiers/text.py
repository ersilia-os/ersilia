import hashlib


class TextIdentifier(object):
    def __init__(self):
        pass

    @staticmethod
    def _is_checksum(text):
        # TODO this method is not secure, for now, it is just a quick check
        if not text.startswith("key"):
            return False
        if len(text) != 32 + 3:
            return False
        text = text[3:]
        if " " in text:
            return False
        if "," in text:
            return False
        if not all(c in "0123456789abcdefABCDEF" for c in text):
            return False
        return True

    def encode(self, text: str) -> str:
        return "key" + hashlib.md5(text.encode("utf-8")).hexdigest()


Identifier = TextIdentifier
