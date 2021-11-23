import os
from ..default import EOS


def tmp_pid_file(model_id):
    """Gets filename to store process ID"""
    return os.path.join(EOS, "tmp", "{0}.pid".format(model_id))
