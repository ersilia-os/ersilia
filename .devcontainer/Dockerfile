FROM python:3.10
ENV PATH="/root/miniconda3/bin:$PATH"
ARG PATH="/root/miniconda3/bin:$PATH"

RUN apt-get update && apt-get upgrade -y

RUN apt-get install -y tzdata
RUN apt-get install -y build-essential wget git curl zip docker

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-py310_23.3.1-0-Linux-x86_64.sh -b \
    && rm -f Miniconda3-py310_23.3.1-0-Linux-x86_64.sh
RUN conda update -y conda && conda init && conda --version

ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache

RUN eval "$(/root/miniconda3/bin/conda shell.bash hook)" &&\
    conda create -y -n ersilia python=3.10 && \
    conda activate ersilia

# Install Ersilia
RUN git clone https://github.com/ersilia-os/ersilia.git
RUN python3 -m pip install -e ersilia/

ENV PATH "$PATH:/root/miniconda3/envs/ersilia/bin"
