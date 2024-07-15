FROM ersiliaos/base-v2:latest AS build
ARG MODEL=eos_identifier
ENV MODEL=$MODEL
RUN ersilia -v fetch $MODEL --from_github
# Install conda pack so we can create a standalone environment from
# the model's conda environment
RUN conda install -c conda-forge -y conda-pack
RUN conda-pack -n $MODEL -o /tmp/env.tar && \
    mkdir /$MODEL && cd /$MODEL && tar xf /tmp/env.tar && \
    rm /tmp/env.tar
RUN /$MODEL/bin/conda-unpack

# Now we can create a new image that only contains
# the ersilia environment, the model environment, 
# and the model itself (as a bentoml bundle)

FROM python:3.7-slim-buster
WORKDIR /root
ARG MODEL=eos_identifier
ENV MODEL=$MODEL

# The following lines ensure that ersilia environment is directly in the path
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PIP_DISABLE_PIP_VERSION_CHECK=1
ENV PATH="/venv/bin:$PATH"

# We install nginx here directly instead of the base image
RUN apt-get update && apt-get install nginx -y

# Copy the model env and ersilia env from the build image
COPY --from=build /$MODEL /$MODEL
COPY --from=build /venv /venv

# Retain the bundled model bento from the build stage
COPY --from=build /root/eos /root/eos
# Copy bentoml artifacts so it doesn't complain about model bento not being found
COPY --from=build /root/bentoml /root/bentoml

COPY --from=build /root/docker-entrypoint.sh docker-entrypoint.sh
COPY --from=build /root/nginx.conf /etc/nginx/sites-available/default

# Patch the python path so we call the correct python binary
RUN /bin/bash <<EOF
set -eux
cd /root/bentoml/repository/$MODEL/*/$MODEL/artifacts/
if [ -f framework/run.sh ]; then
  sed -i -E 's/(.usr.bin.conda.envs(.*bin.*python)(.*))/\2 \3/' framework/run.sh
fi
EOF

RUN chmod + docker-entrypoint.sh
EXPOSE 80
ENTRYPOINT [ "sh", "docker-entrypoint.sh"]
