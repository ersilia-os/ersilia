#!/bin/bash
set -ex
if [ -z "${MODEL}" ];
then
    echo "Model name has not been specified"
    exit 1
fi

ARGS+=( "--bundle_path" "/root/bundles/$MODEL" "--port" "80" )
if [ -n "${ROOT_PATH}" ];
then
    ARGS+=( "--root_path" "${ROOT_PATH}" )
fi

ersilia_model_serve "${ARGS[@]}"

echo "Serving model $MODEL..."
