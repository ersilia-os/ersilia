FROM python:3.7-slim-buster

WORKDIR /root
ARG BUILD_MODEL
ENV PATH=$PATH:/usr/bin/conda/bin
COPY download-miniconda.sh /root/download-miniconda.sh
RUN chmod + /root/download-miniconda.sh
RUN bash /root/download-miniconda.sh && \
    apt-get install -y libxrender1 build-essential curl git && \
    pip install gh && \
    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && \
    apt-get install git-lfs && \
    git lfs install && \
    mkdir -p /usr/bin/conda && \
    bash /root/miniconda.sh -b -u -p /usr/bin/conda && \
    rm /root/miniconda.sh && \
    git clone https://github.com/DhanshreeA/ersilia.git
# Add venv in path and create a virtutalenv to install ersilia in, 
# so we can retain this env in a later image
ENV PATH="/venv/bin:$PATH"
RUN python -m venv /venv && \
    pip install --upgrade pip && \
    cd /root/ersilia && \
    pip install .
