FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="cb2c51dd0ad3df620c4914840c5ef6f5570a5ffd8cfd54cec57d2ffef0a76b08"

# Step 1: Retrieve conda environments

RUN mkdir -p /conda-envs/a4fcaaffb623ea8aef412c66280bd623
COPY envs/environment_minimal.yml /conda-envs/a4fcaaffb623ea8aef412c66280bd623/environment.yaml

RUN mkdir -p /conda-envs/ef25c8d726aebbe9e0ee64fee6c3caa9
COPY envs/esm.yml /conda-envs/ef25c8d726aebbe9e0ee64fee6c3caa9/environment.yaml

RUN mkdir -p /conda-envs/7f88b844a05ae487b7bb6530b5e6a90c
COPY envs/mmsplice.yml /conda-envs/7f88b844a05ae487b7bb6530b5e6a90c/environment.yaml

RUN mkdir -p /conda-envs/dfc51ced08aaeb4cbd3dcd509dec0fc5
COPY envs/regulatorySequence.yml /conda-envs/dfc51ced08aaeb4cbd3dcd509dec0fc5/environment.yaml

RUN mkdir -p /conda-envs/89fe1049cc18768b984c476c399b7989
COPY envs/vep.yml /conda-envs/89fe1049cc18768b984c476c399b7989/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/a4fcaaffb623ea8aef412c66280bd623 --file /conda-envs/a4fcaaffb623ea8aef412c66280bd623/environment.yaml && \
    mamba env create --prefix /conda-envs/ef25c8d726aebbe9e0ee64fee6c3caa9 --file /conda-envs/ef25c8d726aebbe9e0ee64fee6c3caa9/environment.yaml && \
    mamba env create --prefix /conda-envs/7f88b844a05ae487b7bb6530b5e6a90c --file /conda-envs/7f88b844a05ae487b7bb6530b5e6a90c/environment.yaml && \
    mamba env create --prefix /conda-envs/dfc51ced08aaeb4cbd3dcd509dec0fc5 --file /conda-envs/dfc51ced08aaeb4cbd3dcd509dec0fc5/environment.yaml && \
    mamba env create --prefix /conda-envs/89fe1049cc18768b984c476c399b7989 --file /conda-envs/89fe1049cc18768b984c476c399b7989/environment.yaml && \
    mamba clean --all -y
