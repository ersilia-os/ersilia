#!/bin/bash

set -eux

# Change into the BentoML model bundle artifact directory
cd ~/bentoml/repository/$MODEL/*/$MODEL/artifacts/framework

# Patch the python path in run.sh here to use the standalone python binary
# in model environment within root directory that we created using conda-pack
sed -i -n 's/\/usr\/bin\/conda\/envs//p' run.sh