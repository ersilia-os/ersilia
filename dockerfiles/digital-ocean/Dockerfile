FROM ubuntu:22.04 as conda
RUN apt-get update && \
    apt install -y build-essential curl wget git ca-certificates gnupg lsb-release && \
    mkdir -p /etc/apt/keyrings && \
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /etc/apt/keyrings/docker.gpg && \
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null && \
    apt-get update && \
    apt-get install -y docker-ce-cli && \
    mkdir -p ~/miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh && \
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3 && \
    rm -rf ~/miniconda3/miniconda.sh
SHELL ["/bin/bash", "-c"]
ENV PATH=$PATH:/root/miniconda3/bin
COPY setup-ersilia /usr/bin/setup-ersilia
COPY docker-entrypoint.sh /
RUN chmod +x /usr/bin/setup-ersilia
RUN /usr/bin/setup-ersilia
WORKDIR /root
CMD bash
