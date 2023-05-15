apt-get update
apt-get install -y wget
ARCH=$(uname -m)
if [ "$ARCH" = "x86_64" ]; then
    echo "AMD64 architecture detected" &&
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
    echo "$ARCH" > arch.sh
elif [ "$ARCH" = "arm64" ]; then
    echo "ARM64 architecture detected" &&
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O miniconda.sh;
    echo "$ARCH" > arch.sh
elif [ "$ARCH" = "aarch64" ]; then
    echo "AARCH64 architecture detected" &&
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O miniconda.sh;
    echo "$ARCH" > arch.sh
else
    echo "Unsupported architecture: $ARCH";
    echo "$ARCH" > arch.sh
fi;
