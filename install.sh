#!/bin/bash

set -e

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

cd $BASEDIR

# check whether conda and snakemake are available

if [ "$(type conda)" == '' ]
then
    echo 'Conda seems not to be available. Are you sure conda is installed and available in the current $PATH ?';
    exit 1;
fi
#if [ "$(type snakemake)" == '' ]
#then
#    echo 'Snakemake seems not to be available. Are you sure snakemake is installed and available in the current $PATH ?';
#    exit 1;
#fi

# ask which parts of CADD the user wants to install
read -p "Do you want to install the virtual environment with all CADD dependencies via conda? (y)/n " CHOICE
case "$CHOICE" in
    y|Y ) ENV=true;;
    n|N ) ENV=false;;
    * ) ENV=true; echo "Assuming Yes.";;
esac

read -p "Do you want to install CADD v1.4 for GRCh37/hg19? (y)/n " CHOICE
case "$CHOICE" in
    y|Y ) GRCh37=true;;
    n|N ) GRCh37=false;;
    * ) GRCh37=true; echo "Assuming Yes.";;
esac

read -p "Do you want to install CADD v1.4 for GRCh38/hg38? (y)/n " CHOICE
case "$CHOICE" in
    y|Y ) GRCh38=true;;
    n|N ) GRCh38=false;;
    * ) GRCh38=true; echo "Assuming Yes.";;
esac

if [ "$GRCh37" = false ] && [ "$GRCh38" = false ]
then
    echo "You have choosen to not install any of the two available genome builds. Discontinuing installation.";
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
fi

### FILE CONFIGURATION

ANNOTATION_GRCh37="http://krishna.gs.washington.edu/download/CADD/v1.4/annotationsGRCh37.tar.gz"
ANNOTATION_GRCh38="http://krishna.gs.washington.edu/download/CADD/v1.4/annotationsGRCh38.tar.gz"
PRESCORE_GRCh37="http://krishna.gs.washington.edu/download/CADD/v1.4/whole_genome_SNVs_GRCh37.tsv.gz"
PRESCORE_GRCh38="http://krishna.gs.washington.edu/download/CADD/v1.4/whole_genome_SNVs_GRCh38.tsv.gz"
PRESCORE_INCANNO_GRCh37="http://krishna.gs.washington.edu/download/CADD/v1.4/whole_genome_SNVs_inclAnno_GRCh37.tsv.gz"
PRESCORE_INCANNO_GRCh38="http://krishna.gs.washington.edu/download/CADD/v1.4/whole_genome_SNVs_inclAnno_GRCh38.tsv.gz"

### OVERVIEW SELECTION

echo "The following will be loaded: (disk space occupied)"

if [ "$ENV" = true ]
then
    echo " - Setup of the virtual environment including all dependencies (3 GB)."
fi

if [ "$GRCh37" = true ]
then
    if [ "$ANNOTATIONS" = true ]
    then
        echo " - Download CADD annotations for GRCh37 (99 GB)"
    fi

    if [ "$PRESCORE" = true ]
    then
        if [ "INCANNO" = true ]
        then
            echo " - Download prescored SNV inclusive annotations for GRCh37 (500 GB)"
        fi
        if [ "NOANNO" = true ]
        then
            echo " - Download prescored SNV (without annotations) for GRCh37 (80 GB)"
        fi
    fi
fi

if [ "$GRCh38" = true ]
then
    if [ "$ANNOTATIONS" = true ]
    then
        echo " - Download CADD annotations for GRCh38 (194 GB)"
    fi

    if [ "$PRESCORE" = true ]
    then
        if [ "INCANNO" = true ]
        then
            echo " - Download prescored SNV inclusive annotations for GRCh38 (500 GB)"
        fi
        if [ "NOANNO" = true ]
        then
            echo " - Download prescored SNV (without annotations) for GRCh38 (80 GB)"
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

if [ "$ENV" = true ]
then
    echo "Setting up virtual environment"
    conda env create -f src/environment.yml
fi

if [ "$GRCh37" = true ]
then
    if [ "$ANNOTATIONS" = true ]
    then
        echo "Downloading CADD annotations for GRCh37 (99 GB)"
        cd data/annotation/
        wget -c ${ANNOTATION_GRCh37} -O annotationsGRCh37.tar.gz
        wget -c ${ANNOTATION_GRCh37}.md5 -O annotationsGRCh37.tar.gz.md5
        md5sum -c annotationsGRCh37.tar.gz.md5
        echo "Unpacking CADD annotations for GRCh37"
        tar -zxf annotationsGRCh37.tar.gz
        rm annotationsGRCh37.tar.gz
        rm annotationsGRCh37.tar.gz.md5
        cd $OLDPWD
    fi

    if [ "$PRESCORE" = true ]
    then
        if [ "$NOANNO" = true ]
        then
            echo "Downloading prescored SNV without annotations for GRCh37 (80 GB)"
            cd data/prescored/GRCh37/no_anno/
            wget -c ${PRESCORE_GRCh37}
            wget -c ${PRESCORE_GRCh37}.tbi
            wget -c ${PRESCORE_GRCh37}.md5
            md5sum -c *.md5
            cd $OLDPWD
        fi

        if [ "$INCANNO" = true ]
        then
            echo "Downloading prescored SNV inclusive annotations for GRCh37 (500 GB)"
            cd data/prescored/GRCh37/no_anno/
            wget -c ${PRESCORE_INCANNO_GRCh37}
            wget -c ${PRESCORE_INCANNO_GRCh37}.tbi
            wget -c ${PRESCORE_INCANNO_GRCh37}.md5
            md5sum -c *.md5
            cd $OLDPWD
        fi
    fi
fi

if [ "$GRCh38" = true ]
then

    if [ "$ANNOTATIONS" = true ]
    then
        echo "Downloading CADD annotations for GRCh38 (194 GB)"
        cd data/annotation/
        wget -c $ANNOTATION_GRCh38 -O annotationsGRCh38.tar.gz
        wget -c $ANNOTATION_GRCh38.md5 -O annotationsGRCh38.tar.gz.md5
        md5sum -c annotationsGRCh38.tar.gz.md5
        echo "Unpacking CADD annotations for GRCh38"
        tar -zxf annotationsGRCh38.tar.gz
        rm annotationsGRCh38.tar.gz
        rm annotationsGRCh38.tar.gz.md5
        cd $OLDPWD
    fi

    if [ "$PRESCORE" = true ]
    then
        if [ "$NOANNO" = true ]
        then
            echo "Downloading prescored SNV without annotations for GRCh38 (80 GB)"
            cd data/prescored/GRCh38/incl_anno/
            wget -c ${PRESCORE_GRCh38}
            wget -c ${PRESCORE_GRCh38}.tbi
            wget -c ${PRESCORE_GRCh38}.md5
            md5sum -c *.md5
            cd $OLDPWD
        fi

        if [ "$INCANNO" = true ]
        then
            echo "Downloading prescored SNV inclusive annotations for GRCh38 (500 GB)"
            cd data/prescored/GRCh38/incl_anno/
            wget -c ${PRESCORE_INCANNO_GRCh38}
            wget -c ${PRESCORE_INCANNO_GRCh38}.tbi
            wget -c ${PRESCORE_INCANNO_GRCh38}.md5
            md5sum -c *.md5
            cd $OLDPWD
        fi
    fi
fi
