import os
import sys

CHECKSUM_NCHAR = 8
CHECKSUM_FILE = ".conda_checksum"


def read_checksum(path):
    fn = os.path.join(path, CHECKSUM_FILE)
    if not os.path.exists(fn):
        return None
    with open(fn, "r") as f:
        checksum = f.read().strip()
    return checksum


if __name__ == "__main__":
    args = sys.argv
    if len(args) < 2:
        sys.exit()
    model_id = args[1]
    root = os.path.dirname(os.path.realpath(__file__))
    path = os.path.join(root, "dest", model_id)
    checksum = read_checksum(path)
    if checksum is None:
        print(model_id)
    else:
        print(checksum)
