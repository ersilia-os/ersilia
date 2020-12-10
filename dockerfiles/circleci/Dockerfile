# stored in dockerhub as ersiliaos/conda:3.7
FROM ubuntu:20.04
MAINTAINER ersilia

# update apt and get miniconda
RUN apt-get update \
    && apt-get install -y wget \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && apt-get install -y openssh-client git \
    && apt-get install -y gcc

# install miniconda
ENV PATH="/root/miniconda3/bin:$PATH"
RUN mkdir /root/.conda && bash Miniconda3-latest-Linux-x86_64.sh -b

# create conda environment
RUN conda init bash \
    && . ~/.bashrc \
    && conda create --name ersilia python=3.7 \
    && conda activate ersilia \
