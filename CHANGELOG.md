# CHANGELOG

## v1.7.3

- fix that long MNVs can cause sequences of != 500 in regseq annotaton rule and cause a workflow failure ([#90](https://github.com/kircherlab/CADD-scripts/pull/90), [#89](https://github.com/kircherlab/CADD-scripts/issues/89)).

## v1.7.2

- only snakemake >= 8.25.2 supported
- using only conda-forge and bioconda channels (no default anymore)
- new container docker://visze/cadd-scripts-v1_7:0.1.1
- only conda >24.7.1 is allowed (no mamba support anymore)
- VCF2vepVCF.py script fix to extend header. Otherwise regseq will fail using the vcf library
- readme update


## v1.7.1

- containerization
- updating dependencies to resolve conda build issues
- snakemake 8 compatibility
