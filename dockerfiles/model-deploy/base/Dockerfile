FROM python:3.7-slim-buster
WORKDIR /root
ARG BUILD_MODEL
ENV PATH=$PATH:/usr/bin/conda/bin
RUN apt-get update && \
    apt-get install -y build-essential wget curl git nginx && \
    pip install gh && \
    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && \
    apt-get install git-lfs && \
    git lfs install && \
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /root/miniconda.sh && \
    mkdir -p /usr/bin/conda && \
    bash /root/miniconda.sh -b -u -p /usr/bin/conda && \
    rm /root/miniconda.sh && \
    git clone https://github.com/ersilia-os/ersilia.git && \
    cd /root/ersilia && pip install -e .
COPY docker-entrypoint.sh /root/docker-entrypoint.sh
COPY nginx.conf /etc/nginx/sites-available/default
RUN chmod + /root/docker-entrypoint.sh
EXPOSE 80
ENTRYPOINT [ "sh", "/root/docker-entrypoint.sh"]
