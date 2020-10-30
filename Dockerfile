FROM bentoml/model-server
MAINTAINER ersilia

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

WORKDIR /usr/src/

COPY . .

RUN pip install joblib
RUN conda install -c conda-forge rdkit
RUN conda install -c conda-forge biopython
RUN pip install .