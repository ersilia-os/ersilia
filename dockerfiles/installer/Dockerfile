FROM ubuntu:20.04
SHELL ["/bin/bash", "-c"]
WORKDIR /root
RUN apt-get update
RUN apt-get install -y tzdata
RUN apt-get install -y build-essential wget git curl zip docker
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh
COPY download-miniconda.sh /root/download-miniconda.sh
RUN chmod +x /root/download-miniconda.sh
RUN bash /root/download-miniconda.sh
RUN bash /root/miniconda.sh -b -u
RUN eval "$(/root/miniconda3/bin/conda shell.bash hook)" &&\
conda create -y -n ersilia python=3.10 &&\
conda activate ersilia
ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN git clone https://github.com/ersilia-os/ersilia.git
RUN eval "$(/root/miniconda3/bin/conda shell.bash hook)" &&\
    conda activate ersilia &&\
    python -m pip install -e ersilia/. &&\
    conda install -c conda-forge git-lfs &&\
    conda install -c conda-forge gh
RUN rm arch.sh &&\
    rm Miniconda* &&\
    rm download-miniconda.sh &&\
    rm miniconda.sh
ENV PATH "$PATH:$PATH:/root/miniconda3/envs/ersilia/bin"
RUN mv /root/miniconda3/envs/ersilia/bin/ersilia /root/miniconda3/envs/ersilia/bin/ersilia_exec
COPY ./required_files/ersilia /root/miniconda3/envs/ersilia/bin/
RUN chmod +x /root/miniconda3/envs/ersilia/bin/ersilia
ENTRYPOINT ["tail", "-f", "/dev/null"]
