FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="ba86157bf6dc1764a631dfe57d37534afedc875d0531d3ec01927ba2e1b805b4"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: envs/environment_minimal.yml
#   prefix: /conda-envs/919fede82b5b91aac264bb549dd8947a
#   ---
#   name: cadd-minimal
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - ensembl-vep=95.1
#     - python=2.7
#     - scikit-learn=0.20.4
#     - htslib=1.9
#     - pysam=0.11.2
#     - pyvcf=0.6.8
#     - pandas
RUN mkdir -p /conda-envs/919fede82b5b91aac264bb549dd8947a
COPY envs/environment_minimal.yml /conda-envs/919fede82b5b91aac264bb549dd8947a/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/919fede82b5b91aac264bb549dd8947a --file /conda-envs/919fede82b5b91aac264bb549dd8947a/environment.yaml && \
    conda clean --all -y
