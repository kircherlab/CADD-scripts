#!/bin/bash

# CADD_GRCh37-v1.4.sh <input vcf or vcf.gz file>

set -e

FILENAME=${1/*\/}
NAME=${FILENAME/\.*/}
FILEDIR=${1/$FILENAME*/}
FILEFORMAT=${FILENAME/$NAME\./}

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

if [ $FILEFORMAT != "vcf" ] && [ $FILEFORMAT != "vcf.gz" ]
then
    echo "Unknown file format $FILEFORMAT. Make sure you provide a *.vcf or *.vcf.gz file."
    exit 1
fi

cp $1 $BASEDIR/input/GRCh37/$FILENAME

cd $BASEDIR

snakemake --configfile config/CADD1.4-GRCh37.yml --use-conda output/GRCh37/${NAME}.tsv.gz

cd $OLDPWD

mv $BASEDIR/output/GRCh37/${NAME}.tsv.gz ${FILEDIR}${NAME}.tsv.gz
rm $BASEDIR/input/GRCh37/$FILENAME

echo -e "\nScored CADD variants written to file ${FILEDIR}${NAME}.tsv.gz"

exit 0
