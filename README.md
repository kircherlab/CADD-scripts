## Combined Annotation Dependent Depletion (CADD)

CADD is a tool for scoring the deleteriousness of single nucleotide variants as well as insertion/deletions variants in the human genome (currently supported builds: GRCh37/hg19 and GRCh38/hg38).

Details about CADD, including features in the latest version, the different genome builds and how we envision the use case of CADD are described in our latest manuscript:
<blockquote>
    Schubach M, Maass T, Nazaretyan L, Roner S, Kircher M. <br>
    <i>CADD v1.7: Using protein language models, regulatory CNNs and other nucleotide-level scores to improve genome-wide variant predictions.</i><br>
    Nucleic Acids Res. 2023 Nov. doi: <a target="_blank"
        href="https://doi.org/10.1093/nar/gkad989">10.1093/nar/gkad989</a>.<br>
    <!-- PubMed PMID: <a target="_blank" href="http://www.ncbi.nlm.nih.gov/pubmed/TODO">TODO</a>. -->
</blockquote>

And and in our previous publications:
<blockquote>
    Rentzsch P, Schubach M, Shendure J, Kircher M. <br>
    <i>CADD-Splice—improving genome-wide variant effect prediction using deep learning-derived splice scores.</i><br>
    Genome Med. 2021 Feb 22. doi: <a target="_blank" href="https://doi.org/10.1186/s13073-021-00835-9">10.1186/s13073-021-00835-9</a>.<br>
    PubMed PMID: <a target="_blank" href="http://www.ncbi.nlm.nih.gov/pubmed/33618777">33618777</a>.
</blockquote>

<blockquote>
    Rentzsch P, Witten D, Cooper GM, Shendure J, Kircher M. <br>
    <i>CADD: predicting the deleteriousness of variants throughout the human genome.</i><br>
    Nucleic Acids Res. 2018 Oct 29. doi: <a target="_blank" href="http://dx.doi.org/10.1093/nar/gky1016">10.1093/nar/gky1016</a>.<br>
    PubMed PMID: <a target="_blank" href="http://www.ncbi.nlm.nih.gov/pubmed/30371827">30371827</a>.
</blockquote>

The original manuscript describing the method and its features was published by Nature Genetics in 2014:
<blockquote>
    Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J. <br>
    <i>A general framework for estimating the relative pathogenicity of human genetic variants.</i><br>
    Nat Genet. 2014 Feb 2. doi: <a target="_blank" href="http://dx.doi.org/10.1038/ng.2892">10.1038/ng.2892</a>.<br>
    PubMed PMID: <a target="_blank" href="http://www.ncbi.nlm.nih.gov/pubmed/24487276">24487276</a>.
</blockquote>

We provide pre-computed CADD-based scores (C-scores) for all 9 billion possible single nucleotide variants (SNVs) of the reference genome, as well as
all SNV and insertions/deletions variants (InDels) from population-wide whole genome variant releases and enable scoring of short InDels on our website.

Please check our [website for updates and further information](https://cadd.bihealth.org/)

## Offline Installation

This section describes how users can setup CADD version 1.7 on their own system. Please note that this requires between 100 GB - 1 TB of disc space and at least 12 GB of RAM.

### Prerequisite

- conda or mamba

  We recommend to install conda/mamba via [miniforge](https://github.com/conda-forge/miniforge)
```bash
# For example 
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```
- snakemake 7.X (installed via mamba)
```bash
micromamba install -c conda-forge -c bioconda 'snakemake=7'
```

*Note1: If you are using an existing conda installation, please make sure it is [a version >=4.4.0](https://github.com/conda/conda/issues/3200). Make also sure to use snakemake >= 7.32.3 as some command line parameters are not available in earlier versions. *

*Note2: We are using mamba here. In principle it should also work with conda, in that case add `--conda-frontend conda` to line 216 in install.sh`*

*Note3: The recent snakemake 8 is not supported yet. please use the latest snakemake 7.x version*

### Setup

- load/move the zipped CADD archive to its destination folder and unzip.

```bash
unzip CADD.zip
```

- from here, you can either install everything seperately or run the script `install.sh` to install using a brief installation dialog (see below).

#### Install script

This is the easier way of installing CADD, just run:

```
./install.sh
```

You first state which parts you want to install (the environments as well as at least one genome build including annotation tracks are neccessary for a quick start) and the script should manage loading and unpacking the neccessary files.

#### Manual installation

Running CADD depends on four big building blocks (plus the repository containing this README which we assume you already downloaded):

 - snakemake
 - dependencies
 - genome annotations
 - prescored variants

**Installing dependencies**

As of this version, dependencies have to be installed via conda and snakemake. This is because we are using two different enviroments for python2 and python3.

```bash
snakemake test/input.tsv.gz --use-conda --conda-create-envs-only --conda-prefix envs/conda \
        --configfile config/config_GRCh38_v1.7.yml --snakefile Snakefile -c 1
```

Please note that we installing both conda environments in the CADD subdirectory `envs` via `--conda-prefix envs`. If you do not want this behavior (we do this in order to not install the environments in all active directories you run CADD from), adjust or remove this parameter.

**Installing annotations**

Both version of CADD (for the different genome builds) rely on a big number of genomic annotations. Depending on which genome build you require you can get them from our website (be careful where you put them as these are really big files and have identical filenames) via:

```bash
cd data/annotations
# for GRCh37 / hg19
wget -c https://cadd.bihealth.org/download/CADD/v1.7/GRCh37/annotationsGRCh37_v1.7.tar.gz
# for GRCh38 / hg38
wget -c https://cadd.bihealth.org/download/CADD/v1.7/GRCh38/annotationsGRCh38_v1.7.tar.gz
cd $OLDPWD
```

As those files are about 100 and 200 GB in size, downloads can take long (depending on your internet connection). We recommend to setup the process in the background and using a tool (like `wget -c` mentioned above) that allows you to continue an interrupted download.

To make sure you downloaded the files correctly, we recommend downloading md5 hash files from our website (e.g. `wget https://cadd.bihealth.org/download/CADD/v1.7/GRCh38/annotationsGRCh38_v1.7.tar.gz.md5`) and checking for completeness (via `md5sum -c annotationsGRCh38_v1.7.tar.gz.md5`).

The annotation files are finally put in the folder `data/annotations` and unpacked:

```bash
cd data/annotations
# for GRCh37 / hg19
tar -zxvf annotationsGRCh37_v1.7.tar.gz
# for GRCh38 / hg38
tar -zxvf annotationsGRCh38_v1.7.tar.gz
cd $OLDPWD
```

**Installing prescored files**

At this point you are ready to go, but if you want a faster version of CADD, you can download the prescored files from our website (see section Downloads for a list of available files). Please note that these files can be very big. The files are (together with their respective tabix indices) put in the folders `no_anno` or `incl_anno` depending on the file under `data/prescored/${GENOME_BUILD}_${VERSION}/` and will be automatically detected by the `CADD.sh` script.

### Running CADD

You run CADD via the script `CADD.sh` which technically only requieres an either vcf or vcf.gz input file as last argument. You can further specify the genome build via `-g`, CADD version via `-v` (deprecated since version v1.6), request a fully annotated output (`-a` flag) and specify a seperate output file via `-o` (else inputfile name `.tsv.gz` is used). I.e:

```bash
./CADD.sh test/input.vcf

./CADD.sh -a -g GRCh37 -o output_inclAnno_GRCh37.tsv.gz test/input.vcf
```

You can test whether your CADD is set up properly by comparing to the example files in the `test` directory.

### Update

Version 1.7 includes some additions in features to v1.6 and we refactured the Snakemake workflow of CADD-scripts including config files. The new models for v1.7 are extended by a protein language model, regulatory effect prediction with a CNN, Zonoomia and Aparent2 scores and an update to the lates gencode 110 version. Because of the refactored workflow with additinal config settings you cannot use it with previous CADD versions. If you are still using those versions, please use [this repository for version 1.6](https://github.com/kircherlab/CADD-scripts/archive/v1.6.post1.zip) or [this repository for v1.5 and v1.4](https://github.com/kircherlab/CADD-scripts/archive/CADD1.5.zip).

Version 1.6 includes some changes in comparison to v1.5. Next to the obvious switch of the pipeline into a Snakemake workflow which became necessary due to the ongoin issues with `conda activate`, the new models for v1.6 are extended by more specialized annotations for splicing variants, as well as a few minor changes in some other annotations (most prominent: fixed gerp for GRCh38) and changes in consequence categories which make this scripts incompatible with CADD v1.4 and v1.5. If you are still using those version, please use [version 1.5 of this repository](https://github.com/kircherlab/CADD-scripts/archive/CADD1.5.zip).

```

## Copyright
Copyright (c) University of Washington, Hudson-Alpha Institute for
Biotechnology and Berlin Institute of Health at Charité - Universitätsmedizin
Berlin 2013-2023. All rights reserved.

Permission is hereby granted, to all non-commercial users and licensees of CADD
(Combined Annotation Dependent Framework, licensed by the University of
Washington) to obtain copies of this software and associated documentation
files (the "Software"), to use the Software without restriction, including
rights to use, copy, modify, merge, and distribute copies of the Software. The
above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
