FROM mambaorg/micromamba:1.5.10

# Install conda env
COPY TEforest.yml /tmp/TEforest.yml
RUN micromamba create -y -f /tmp/TEforest.yml && micromamba clean -a -y

ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=/opt/micromamba/envs/TEforest/bin:$PATH
ENV TEFOREST_HOME=/opt/teforest

# Copy TEforest code + models into the image
COPY . /opt/teforest
WORKDIR /opt/teforest

ENTRYPOINT ["python", "/opt/teforest/TEforest.py"]
