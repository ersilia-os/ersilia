#!/bin/bash
set -e

# Determine architecture
ARCH=$(uname -m)

case "$ARCH" in
    x86_64 | amd64)
        echo "AMD64 architecture detected"
        URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
        ;;
    arm64 | aarch64)
        echo "ARM64 architecture detected"
        URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh"
        ;;
    *)
        echo "Unsupported architecture: $ARCH"
        exit 1
        ;;
esac

# Download Miniforge installer
wget "$URL" -O miniconda.sh
echo "$ARCH" > arch.sh