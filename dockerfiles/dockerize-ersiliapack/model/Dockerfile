FROM base-image
ARG MODEL=eos_identifier
ENV MODEL=$MODEL
WORKDIR /root
COPY ./$MODEL /root/$MODEL
RUN mkdir /root/bundles && rm -rf /root/$MODEL/.git && rm -rf /root/$MODEL/.github
RUN ersilia_model_pack --repo_path $MODEL --bundles_repo_path /root 
    