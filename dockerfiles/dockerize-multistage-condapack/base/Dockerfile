# Stage 1: Build environment
FROM python:3.10-slim-buster AS builder

# Set working directory
WORKDIR /root

# Install build dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    curl \
    wget \
    build-essential \
    git \
    git-lfs && \
    git lfs install && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Download and install Miniforge (lightweight alternative to Miniconda)
COPY download-miniconda.sh .
RUN chmod +x download-miniconda.sh && \
    ./download-miniconda.sh && \
    bash miniconda.sh -b -u -p /usr/bin/conda && \
    rm miniconda.sh

# Clone Ersilia repository
RUN git clone https://github.com/ersilia-os/ersilia.git

# Create a virtual environment and install Ersilia
RUN python -m venv /venv && \
    /venv/bin/pip install --no-cache-dir --upgrade pip && \
    /venv/bin/pip install --no-cache-dir /root/ersilia

# Stage 2: Runtime environment
FROM python:3.10-slim-buster

WORKDIR /root

# Copy virtual environment from builder stage
COPY --from=builder /venv /venv

# Copy essential files
COPY docker-entrypoint.sh docker-entrypoint.sh
COPY nginx.conf nginx.conf

ENV PATH="/venv/bin:$PATH"

# Set the entry point
ENTRYPOINT ["./docker-entrypoint.sh"]
