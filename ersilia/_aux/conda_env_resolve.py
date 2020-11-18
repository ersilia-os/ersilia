# imports
import sys
import os
import pip

# get argument
args = sys.argv
if len(args) < 2:
    sys.exit()
model_id = args[1]

# check that environment.yml file exists
env_yml = os.path.join('dest', model_id, 'environment.yml')
if not os.path.exists(env_yml):
    sys.exit()

# install hashlib if necessary
try:
    import hashlib
except:
    pip.main(['install', 'hashlib'])
    import hashlib

# md5 checksum of the config yml file
with open(env_yml, "r") as f:
    R = []
    for r in f:
        R += [r]
    S = R[1:-2]
    text = ''.join(S)
print(text)
checksum = hashlib.md5(text.encode('utf-8')).hexdigest()
with open(env_yml, "w") as f:
    f.write('name: %s\n' % checksum)
    for s in S:
        f.write(s)
    prefix = os.path.join('/'.join(R[-2].split('prefix: ')[1].split('/')[:-1]), checksum)
    f.write('prefix: %s\n' % prefix)
    f.write('\n')

print(checksum)
