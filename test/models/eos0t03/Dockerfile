FROM bentoml/model-server:0.11.0-py37

RUN conda install -c conda-forge rdkit

WORKDIR /repo
COPY . /repo
