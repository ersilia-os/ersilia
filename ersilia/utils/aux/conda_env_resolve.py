import sys
import os
import pip
try:
    import hashlib
except ImportError:
    pip.main(['install', 'hashlib'])
    import hashlib


def checksum_from_conda_yaml_file(env_yml, overwrite):
    with open(env_yml, "r") as f:
        name_idx = None
        pref_idx = None
        R = []
        for i, r in enumerate(f):
            if r[:4] == "name":
                name_idx = i
            if r[:6] == "prefix":
                pref_idx = i
            R += [r]
        S = R[(name_idx+1):pref_idx]
        text = ''.join(S)
    checksum = hashlib.md5(text.encode('utf-8')).hexdigest()
    if overwrite:
        with open(env_yml, "w") as f:
            f.write('name: %s\n' % checksum)
            for s in S:
                f.write(s)
            prefix = os.path.join(os.sep.join(R[pref_idx].split('prefix: ')[1].split(os.sep)[:-1]), checksum)
            f.write('prefix: %s\n' % prefix)
            f.write('\n')
    return checksum


if __name__ == "__main__":
    args = sys.argv
    if len(args) < 2:
        sys.exit()
    model_id = args[1]
    env_yml = os.path.join('dest', model_id, 'environment.yml')
    if not os.path.exists(env_yml):
        sys.exit()
    print(checksum_from_conda_yaml_file(model_id, True))
