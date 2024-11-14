FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="2c2e74f6a7392d16830a1aa1a341e3bd0565e62e3bb0bc4932c8008320fbe806"

# Step 1: Retrieve conda environments

RUN mkdir -p /conda-envs/fe7049a03d989a2faa6d7688f5987d09
COPY envs/environment_minimal.yml /conda-envs/fe7049a03d989a2faa6d7688f5987d09/environment.yaml

RUN mkdir -p /conda-envs/10a710eb7e199e21064fc58a16f9ea4f
COPY envs/esm.yml /conda-envs/10a710eb7e199e21064fc58a16f9ea4f/environment.yaml

RUN mkdir -p /conda-envs/27e6f085c62573d05046d5330f07505d
COPY envs/mmsplice.yml /conda-envs/27e6f085c62573d05046d5330f07505d/environment.yaml

RUN mkdir -p /conda-envs/f675d936396bfaffdd96e0282b787268
COPY envs/regulatorySequence.yml /conda-envs/f675d936396bfaffdd96e0282b787268/environment.yaml

RUN mkdir -p /conda-envs/0b37582dec236e35970e555dc2c0f0b5
COPY envs/vep.yml /conda-envs/0b37582dec236e35970e555dc2c0f0b5/environment.yaml

# Step 2: Generate conda environments

RUN conda env create --no-default-packages --prefix /conda-envs/fe7049a03d989a2faa6d7688f5987d09 --file /conda-envs/fe7049a03d989a2faa6d7688f5987d09/environment.yaml && \
    conda env create --no-default-packages --prefix /conda-envs/10a710eb7e199e21064fc58a16f9ea4f --file /conda-envs/10a710eb7e199e21064fc58a16f9ea4f/environment.yaml && \
    conda env create --no-default-packages --prefix /conda-envs/27e6f085c62573d05046d5330f07505d --file /conda-envs/27e6f085c62573d05046d5330f07505d/environment.yaml && \
    conda env create --no-default-packages --prefix /conda-envs/f675d936396bfaffdd96e0282b787268 --file /conda-envs/f675d936396bfaffdd96e0282b787268/environment.yaml && \
    conda env create --no-default-packages --prefix /conda-envs/0b37582dec236e35970e555dc2c0f0b5 --file /conda-envs/0b37582dec236e35970e555dc2c0f0b5/environment.yaml && \
    conda clean --all -y
