FROM bentoml/model-server:0.11.0-py37

RUN conda install -c conda-forge rdkit=2020.03

WORKDIR /repo
COPY . /repo
