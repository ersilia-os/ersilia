FROM ersiliaos/base-v4:latest
ARG MODEL=eos4zfy
ENV MODEL=$MODEL
RUN pip install --no-cache-dir packaging
RUN ersilia -v fetch $MODEL
