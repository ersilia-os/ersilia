import re
import os


def model_id_from_path(path):
    regex = re.compile(r"eos[0-9][a-z0-9]{3}")
    path = os.path.abspath(path)
    model_ids = sorted(set(regex.findall(path)))
    if len(model_ids) == 1:
        return model_ids[0]
    else:
        return None
