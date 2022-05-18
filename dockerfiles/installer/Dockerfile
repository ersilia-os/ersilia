FROM ubuntu:20.04
RUN apt-get update
RUN apt-get install -y tzdata
RUN apt-get install -y build-essential wget git curl zip docker
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.11.0-Linux-x86_64.sh
RUN chmod +x Miniconda3-py37_4.11.0-Linux-x86_64.sh
RUN bash Miniconda3-py37_4.11.0-Linux-x86_64.sh -b
RUN eval "$(/root/miniconda3/bin/conda shell.bash hook)" &&\
conda create -y -n ersilia python=3.7 &&\
conda activate ersilia
ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN git clone https://github.com/ersilia-os/ersilia.git
RUN eval "$(/root/miniconda3/bin/conda shell.bash hook)" &&\
conda activate ersilia &&\
python -m pip install -e ersilia/.
ENV PATH "$PATH:$PATH:/root/miniconda3/envs/ersilia/bin"
