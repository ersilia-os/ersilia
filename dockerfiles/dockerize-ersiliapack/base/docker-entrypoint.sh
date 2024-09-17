#!/bin/bash
set -ex
if [ -z "${MODEL}" ];
then
    echo "Model name has not been specified"
    exit 1
fi
ersilia_model_serve --bundle_path /root/bundles/$MODEL --port 3000
echo "Serving model $MODEL..."
nginx -g 'daemon off;'
