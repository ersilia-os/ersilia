#!/bin/bash
set -ex
if [ -z "${MODEL}" ];
then
    echo "Model name has not been specified"
    exit 1
fi
ersilia serve -p 3000 $MODEL
echo "Serving model $MODEL..."
nginx -g 'daemon off;'
