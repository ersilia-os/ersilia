# Stage 1: Build environment
FROM ersiliaos/base-v4:latest AS build

# Install git and other necessary packages before running ersilia fetch
RUN apt-get update && \
  apt-get install -y --no-install-recommends git git-lfs && git lfs install  && \
  rm -rf /var/lib/apt/lists/*

# Fetch model and pack environment
ARG MODEL=eos4zfy
ENV MODEL=$MODEL

RUN pip install --no-cache-dir packaging

RUN ersilia -v fetch $MODEL --from_github && \
  conda install -c conda-forge -y conda-pack && \
  conda-pack -n $MODEL -o /tmp/env.tar && \
  mkdir /$MODEL && tar xf /tmp/env.tar -C /$MODEL && \
  rm /tmp/env.tar && \
  /$MODEL/bin/conda-unpack

# Install necessary Python dependencies (packaging)
RUN /$MODEL/bin/pip install --no-cache-dir packaging

# Stage 2: Final runtime image
FROM python:3.10-slim-buster

# Set working directory
WORKDIR /root

# Environment variables
ARG MODEL=eos_identifier
ENV MODEL=$MODEL
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PIP_DISABLE_PIP_VERSION_CHECK=1
ENV PATH="/venv/bin:$PATH"

# Install runtime dependencies
RUN apt-get update && \
  apt-get install -y --no-install-recommends nginx && \
  rm -rf /var/lib/apt/lists/*

# Copy model environment and dependencies from the build stage
COPY --from=build /$MODEL /$MODEL
COPY --from=build /venv /venv

# Copy model artifacts and scripts
COPY --from=build /root/eos /root/eos
COPY --from=build /root/bentoml /root/bentoml
COPY --from=build /root/docker-entrypoint.sh docker-entrypoint.sh
COPY --from=build /root/nginx.conf /etc/nginx/sites-available/default

# Patch Python path for correct binary usage
RUN sed -i -E 's|/usr/bin/conda/envs/[^/]+/bin/python|/venv/bin/python|' \
  /root/bentoml/repository/$MODEL/*/$MODEL/artifacts/framework/run.sh || true

# Make entrypoint script executable
RUN chmod +x docker-entrypoint.sh

# Expose port for service
EXPOSE 80

# Set entrypoint
ENTRYPOINT [ "sh", "docker-entrypoint.sh" ]
