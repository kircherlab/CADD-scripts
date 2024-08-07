#!/bin/bash

usage="$(basename "$0") [-b] -- CADD installer, version 1.7

where:
    -h  show this help text
    -b  use kircherlab.bihealth.org as download server instead of krishna.gs.washington.edu

All further settings will be spcified in a yes/no dialog. 
    "

### SET FILE SERVER
DOWNLOAD_LOCATION="https://krishna.gs.washington.edu/download/CADD"

while getopts ':hb' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    b) DOWNLOAD_LOCATION=https://kircherlab.bihealth.org/download/CADD
       echo "Using kircherlab.bihealth.org as download server"
       ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND-1))

set -e

echo "CADD-v1.7 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health at Charité Universitätsmedizin Berlin 2013-2023. All rights reserved."
echo ""

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export CADD=$(dirname "$SCRIPT")

cd $BASEDIR

# check whether conda and snakemake are available

if [ "$(type conda)" == '' ] && [ "$(type mamba)" == '' ] && [ "$(type micromamba)" == '' ]
then
    echo 'Conda/mamba seems not to be available. Are you sure conda is installed and available in the current $PATH ?';
fi

if [ "$(type snakemake)" == '' ]
then
    echo 'Snakemake seems not to be available. Are you sure snakemake is installed and available in the current $PATH ?';
    exit 1;
fi

echo "The following questions will quide you through selecting the files and dependencies needed for CADD."
echo "After this, you will see an overview of the selected files before the download and installation starts."
echo "Please note, that for successfully running CADD locally, you will need the conda environment and at least one set of annotations."
echo ""

# ask which parts of CADD the user wants to install
read -p "For conda/mamba runs only (running without apptainer): Do you want to install the virtual environments with all CADD dependencies via conda? y/(n) " CHOICE
case "$CHOICE" in
    y|Y ) ENV=true;;
    n|N ) ENV=false;;
    * ) ENV=false; echo "Assuming no.";;
esac

read -p "Do you want to install CADD v1.7 for GRCh37/hg19? (y)/n " CHOICE
case "$CHOICE" in
    y|Y ) GRCh37=true;;
    n|N ) GRCh37=false;;
    * ) GRCh37=true; echo "Assuming Yes.";;
esac

read -p "Do you want to install CADD v1.7 for GRCh38/hg38? (y)/n " CHOICE
case "$CHOICE" in
    y|Y ) GRCh38=true;;
    n|N ) GRCh38=false;;
    * ) GRCh38=true; echo "Assuming Yes.";;
esac

if [ "$GRCh37" = false ] && [ "$GRCh38" = false ]
then
    echo "You have choosen to not install any of the available CADD models. Discontinuing installation.";
    exit 0;
fi

read -p "Do you want to load annotations (Annotations can also be downloaded manually from the website)? (y)/n " CHOICE
case "$CHOICE" in
    y|Y ) ANNOTATIONS=true;;
    n|N ) ANNOTATIONS=false;;
    * ) ANNOTATIONS=true; echo "Assuming Yes.";;
esac

read -p "Do you want to load prescored variants (Makes SNV calling faster. Can also be loaded/installed later.)? y/(n) " CHOICE
case "$CHOICE" in
    y|Y ) PRESCORE=true;;
    n|N ) PRESCORE=false;;
    * ) PRESCORE=false; echo "Assuming No.";;
esac

if [ "$PRESCORE" = true ]
then
    read -p "Do you want to load prescored variants for scoring with annotations (Warning: These files are very big)? y/(n) " CHOICE
    case "$CHOICE" in
        y|Y ) INCANNO=true;;
        n|N ) INCANNO=false;;
        * ) INCANNO=false; echo "Assuming No.";;
    esac
    read -p "Do you want to load prescored variants for scoring without annotations? y/(n) " CHOICE
    case "$CHOICE" in
        y|Y ) NOANNO=true;;
        n|N ) NOCANNO=false;;
        * ) NOANNO=false; echo "Assuming No.";;
    esac
    read -p "Do you also want to load prescored InDels? We provide scores for well known InDels from sources like ClinVar, gnomAD/TOPMed etc. y/(n) " CHOICE
    case "$CHOICE" in
        y|Y ) INDELS=true;;
        n|N ) INDELS=false;;
        * ) INDELS=false; echo "Assuming No.";;
    esac
fi

### FILE CONFIGURATION
ANNOTATION_GRCh37="$DOWNLOAD_LOCATION/v1.7/GRCh37/GRCh37_v1.7.tar.gz"
ANNOTATION_GRCh38="$DOWNLOAD_LOCATION/v1.7/GRCh38/GRCh38_v1.7.tar.gz"
PRESCORE_GRCh37="$DOWNLOAD_LOCATION/v1.7/GRCh37/whole_genome_SNVs.tsv.gz"
PRESCORE_GRCh38="$DOWNLOAD_LOCATION/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
PRESCORE_INCANNO_GRCh37="$DOWNLOAD_LOCATION/v1.7/GRCh37/whole_genome_SNVs_inclAnno.tsv.gz"
PRESCORE_INCANNO_GRCh38="$DOWNLOAD_LOCATION/v1.7/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz"
PRESCORE_GRCh37_INDEL="$DOWNLOAD_LOCATION/v1.7/GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz"
PRESCORE_GRCh38_INDEL="$DOWNLOAD_LOCATION/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz"
PRESCORE_INCANNO_GRCh37_INDEL="$DOWNLOAD_LOCATION/v1.7/GRCh37/gnomad.genomes-exomes.r4.0.indel_inclAnno.tsv.gz"
PRESCORE_INCANNO_GRCh38_INDEL="$DOWNLOAD_LOCATION/v1.7/GRCh38/gnomad.genomes.r4.0.indel_inclAnno.tsv.gz"

### OVERVIEW SELECTION

echo ""
echo "The following will be loaded: (disk space occupied)"

if [ "$ENV" = true ]
then
    echo " - Setup of the virtual environments including all dependencies for CADD v1.7 (16 GB)."
fi

if [ "$GRCh37" = true ]
then
    if [ "$ANNOTATIONS" = true ]
    then
        echo " - Download CADD annotations for GRCh37-v1.7 (261 GB)"
    fi

    if [ "$PRESCORE" = true ]
    then
        if [ "$INCANNO" = true ]
        then
            echo " - Download prescored SNV inclusive annotations for GRCh37-v1.7 (562 GB)"
            if [ "$INDELS" = true ]
            then
                echo " - Download prescored InDels inclusive annotations for GRCh37-v1.7 (10 GB)"
            fi
        fi
        if [ "$NOANNO" = true ]
        then
            echo " - Download prescored SNV (without annotations) for GRCh37-v1.7 (79 GB)"
            if [ "$INDELS" = true ]
            then
                echo " - Download prescored InDels (without annotations) for GRCh37-v1.7 (1.4 GB)"
            fi
        fi
    fi
fi

if [ "$GRCh38" = true ]
then
    if [ "$ANNOTATIONS" = true ]
    then
        echo " - Download CADD annotations for GRCh38-v1.7 (336 GB)"
    fi

    if [ "$PRESCORE" = true ]
    then
        if [ "$INCANNO" = true ]
        then
            echo " - Download prescored SNV inclusive annotations for GRCh38-v1.7 (625 GB)"
            if [ "$INDELS" = true ]
            then
                echo " - Download prescored InDels inclusive annotations for GRCh38-v1.7 (11 GB)"
            fi
        fi
        if [ "$NOANNO" = true ]
        then
            echo " - Download prescored SNV (without annotations) for GRCh38-v1.7 (81 GB)"
            if [ "$INDELS" = true ]
            then
                echo " - Download prescored InDels (without annotations) for GRCh38-v1.7 (1.2 GB)"
            fi
        fi
    fi
fi

echo "Please make sure you have enough disk space available."

read -p "Ready to continue? (y)/n " CHOICE
case "$CHOICE" in
    n|N ) echo "You canceled the installation."; exit 0;;
    * ) echo "Starting installation. This will take some time.";;
esac

### INSTALLATION

# download a file and it index and check both md5 sums
function download_variantfile()
{
    echo $1
    wget -c $2
    wget -c $2.tbi
    wget $2.md5
    wget $2.tbi.md5
    md5sum -c *.md5
    rm *.md5
}

if [ "$GRCh37" = true ]
then
    if [ "$ANNOTATIONS" = true ]
    then
        echo "Downloading CADD annotations for GRCh37-v1.7 (261 GB)"
        mkdir -p data/annotations/
        cd data/annotations/
        wget -c ${ANNOTATION_GRCh37} -O GRCh37_v1.7.tar.gz
        wget ${ANNOTATION_GRCh37}.md5 -O GRCh37_v1.7.tar.gz.md5
        md5sum -c GRCh37_v1.7.tar.gz.md5
        echo "Unpacking CADD annotations for GRCh37-v1.7"
        tar -zxf GRCh37_v1.7.tar.gz
        rm GRCh37_v1.7.tar.gz
        rm GRCh37_v1.7.tar.gz.md5
        cd $OLDPWD
    fi

    if [ "$PRESCORE" = true ]
    then
        if [ "$NOANNO" = true ]
        then
            mkdir -p data/prescored/GRCh37_v1.7/no_anno/
            cd data/prescored/GRCh37_v1.7/no_anno/
            download_variantfile "Downloading prescored SNV without annotations for GRCh37-v1.7 (78 GB)" ${PRESCORE_GRCh37}
            if [ "$INDELS" = true ]
            then
                download_variantfile "Downloading prescored InDels without annotations for GRCh37-v1.7 (1.4 GB)" ${PRESCORE_GRCh37_INDEL}
            fi
            cd $OLDPWD
        fi

        if [ "$INCANNO" = true ]
        then
            mkdir -p data/prescored/GRCh37_v1.7/incl_anno/
            cd data/prescored/GRCh37_v1.7/incl_anno/
            download_variantfile "Downloading prescored SNV inclusive annotations for GRCh37-v1.7 (558 GB)" ${PRESCORE_INCANNO_GRCh37}
            if [ "$INDELS" = true ]
            then
                download_variantfile "Downloading prescored InDels inclusive annotations for GRCh37-v1.7 (11 GB)" ${PRESCORE_INCANNO_GRCh37_INDEL}
            fi
            cd $OLDPWD
        fi
    fi
fi

if [ "$GRCh38" = true ]
then

    if [ "$ANNOTATIONS" = true ]
    then
        echo "Downloading CADD annotations for GRCh38-v1.7 (336 GB)"
        mkdir -p data/annotations/
        cd data/annotations/
        wget -c $ANNOTATION_GRCh38 -O GRCh38_v1.7.tar.gz
        wget $ANNOTATION_GRCh38.md5 -O GRCh38_v1.7.tar.gz.md5
        md5sum -c GRCh38_v1.7.tar.gz.md5
        echo "Unpacking CADD annotations for GRCh38-v1.7"
        tar -zxf GRCh38_v1.7.tar.gz
        rm GRCh38_v1.7.tar.gz
        rm GRCh38_v1.7.tar.gz.md5
        cd $OLDPWD
    fi

    if [ "$PRESCORE" = true ]
    then
        if [ "$NOANNO" = true ]
        then
            mkdir -p data/prescored/GRCh38_v1.7/no_anno/
            cd data/prescored/GRCh38_v1.7/no_anno/
            download_variantfile "Downloading prescored SNV without annotations for GRCh38-v1.7 (81 GB)" ${PRESCORE_GRCh38}
            if [ "$INDELS" = true ]
            then
                download_variantfile "Downloading prescored InDels without annotations for GRCh38-v1.7 (1.2 GB)" ${PRESCORE_GRCh38_INDEL}
            fi
            cd $OLDPWD
        fi

        if [ "$INCANNO" = true ]
        then
            mkdir -p data/prescored/GRCh38_v1.7/incl_anno/
            cd data/prescored/GRCh38_v1.7/incl_anno/
            download_variantfile "Downloading prescored SNV inclusive annotations for GRCh38-v1.7 (621 GB)" ${PRESCORE_INCANNO_GRCh38}
            if [ "$INDELS" = true ]
            then
                download_variantfile "Downloading prescored InDels inclusive annotations for GRCh38-v1.7 (12 GB)" ${PRESCORE_INCANNO_GRCh38_INDEL}
            fi
            cd $OLDPWD
        fi
    fi
fi

if [ "$ENV" = true ]
then
    echo "Setting up virtual environments for CADD v1.7"
    snakemake test/input.tsv.gz --use-conda --conda-create-envs-only --conda-prefix envs/conda \
        --cores 1 --configfile config/config_GRCh38_v1.7.yml --snakefile Snakefile
fi
