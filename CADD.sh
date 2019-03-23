#!/bin/bash

usage="$(basename "$0") [-o <outfile>] [-g <genomebuild>] [-v <caddversion>] [-a] <infile>  -- CADD version 1.5

where:
    -h  show this help text
    -o  out tsv.gz file (generated from input file name if not set)
    -g  genome build (supported are GRCh37 and GRCh38 [default: GRCh38])
    -v  CADD version (either v1.4 or v1.5 [default: v1.5])
    -a  include annotation in output
        input vcf of vcf.gz file (required)"

unset OPTARG
unset OPTIND

GENOMEBUILD="GRCh38"
ANNOTATION=false
OUTFILE=""
VERSION="v1.5"
while getopts ':ho:g:v:a' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    o) OUTFILE=$OPTARG
       ;;
    g) GENOMEBUILD=$OPTARG
       ;;
    v) VERSION=$OPTARG
       ;;
    a) ANNOTATION=true
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND-1))

INFILE=$1

echo "CADD-v1.5 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health 2013-2019. All rights reserved."

set -ueo pipefail

### Configuring all the paths

FILENAME=$(basename $INFILE)
NAME=${FILENAME%\.vcf*}
FILEDIR=$(dirname $INFILE)
FILEFORMAT=${FILENAME#$NAME\.}

SCRIPT=$(readlink -f "$0")
export CADD=$(dirname "$SCRIPT")

if [ "$FILEFORMAT" != "vcf" ] && [ "$FILEFORMAT" != "vcf.gz" ]
then
    echo "Unknown file format $FILEFORMAT. Make sure you provide a *.vcf or *.vcf.gz file."
    exit 1
fi

if [ "$OUTFILE" == "" ]
then
    OUTFILE=$FILEDIR/$NAME.tsv.gz
fi

if [ "$GENOMEBUILD" != "GRCh38" ] && [ "$GENOMEBUILD" != "GRCh37" ]
then
    echo "Unknown/Unsupported genome build $GENOMEBUILD. CADD currently only supports GRCh37 and GRCh38."
    exit 1
fi

if [ "$VERSION" != "v1.4" ] && [ "$VERSION" != "v1.5" ]
then
    echo "Unknown/Unsupported CADD version $VERSION. This script currently only supports v1.4 and v1.5."
    exit 1
fi

if [ "$VERSION" == "v1.5" ] && [ "$GENOMEBUILD" == "GRCh37" ]
then
    echo "Please note that CADD scores for GRCh37 version v1.5 are the same as in v1.4."
    VERSION="v1.4"
fi

if [ "$ANNOTATION" = 'true' ]
then
    ANNO_FOLDER="incl_anno"
else
    ANNO_FOLDER="no_anno"
fi

# Pipeline configuration
PRESCORED_FOLDER=$CADD/data/prescored/${GENOMEBUILD}_${VERSION}/$ANNO_FOLDER
REFERENCE_CONFIG=$CADD/config/references_${GENOMEBUILD}_${VERSION}.cfg
IMPUTE_CONFIG=$CADD/config/impute_${GENOMEBUILD}_${VERSION}.cfg
MODEL=$CADD/data/models/$GENOMEBUILD/CADD${VERSION}-$GENOMEBUILD.mod
CONVERSION_TABLE=$CADD/data/models/$GENOMEBUILD/conversionTable_CADD${VERSION}-$GENOMEBUILD.txt

# determine VEP database version
DBVERSION=92
if [ "$GENOMEBUILD" == "GRCh38" ] && [ "$VERSION" == "v1.5" ]
then
    DBVERSION=95
fi

# Setup temporary folder that is removed reliably on exit and is outside of
# the CADD-scripts directory.
TMP_FOLDER=$(mktemp -d)
trap "rm -rf $TMP_FOLDER" ERR EXIT

# Temp files
TMP_PRE=$TMP_FOLDER/$NAME.pre.tsv.gz
TMP_VCF=$TMP_FOLDER/$NAME.vcf
TMP_ANNO=$TMP_FOLDER/$NAME.anno.tsv.gz
TMP_IMP=$TMP_FOLDER/$NAME.csv.gz
TMP_NOV=$TMP_FOLDER/$NAME.nov.tsv.gz

mkdir -p $TMP_FOLDER

### Pipeline

# Loading the environment
if [ "$VERSION" == "v1.4" ]
then
    source activate cadd-env
else
    source activate cadd-env-v1.5
fi

# File preparation
if [ "$FILEFORMAT" == "vcf" ]
then
    cat $INFILE \
    | python $CADD/src/scripts/VCF2vepVCF.py \
    | sort -k1,1 -k2,2n -k3,3 -k4,4 \
    | uniq > $TMP_VCF 
else
    zcat $INFILE \
    | python $CADD/src/scripts/VCF2vepVCF.py \
    | sort -k1,1 -k2,2n -k3,3 -k4,4 \
    | uniq > $TMP_VCF
fi

# Prescoring
echo '## Prescored variant file' | gzip -c > $TMP_PRE;
if [ -d $PRESCORED_FOLDER ]
then
    for PRESCORED in $(ls $PRESCORED_FOLDER/*.tsv.gz)
    do
        cat $TMP_VCF \
        | python $CADD/src/scripts/extract_scored.py --header \
            -p $PRESCORED --found_out=$TMP_PRE.tmp \
        > $TMP_VCF.tmp;
        gzip -c $TMP_PRE.tmp >> $TMP_PRE
        mv $TMP_VCF.tmp $TMP_VCF;
    done;
    rm $TMP_PRE.tmp
fi

# Variant annotation
cat $TMP_VCF \
| vep --quiet --cache --buffer 1000 --no_stats --offline --vcf \
    --dir $CADD/data/annotations/${GENOMEBUILD}_${VERSION}/vep \
    --species homo_sapiens --db_version=$DBVERSION \
    --assembly $GENOMEBUILD --regulatory --sift b \
    --polyphen b --per_gene --ccds --domains --numbers --canonical \
    --total_length --force_overwrite --format vcf --output_file STDOUT \
    --warning_file STDERR \
| python $CADD/src/scripts/annotateVEPvcf.py -c $REFERENCE_CONFIG \
| gzip -c > $TMP_ANNO
rm $TMP_VCF

# Imputation
zcat $TMP_ANNO \
| python $CADD/src/scripts/trackTransformation.py -b \
            -c $IMPUTE_CONFIG -o $TMP_IMP --noheader;

# Score prediction
python $CADD/src/scripts/predictSKmodel.py \
    -i $TMP_IMP -m $MODEL -a $TMP_ANNO \
| python $CADD/src/scripts/max_line_hierarchy.py --all \
| python $CADD/src/scripts/appendPHREDscore.py -t $CONVERSION_TABLE \
| gzip -c > $TMP_NOV;
rm $TMP_ANNO
rm $TMP_IMP

if [ "$ANNOTATION" = 'false' ]
then
    if [ "$GENOMEBUILD" == "GRCh38" ]
    then
        COLUMNS="1-4,124,125"
    else
        COLUMNS="1-4,106,107"
    fi
    zcat $TMP_NOV | cut -f $COLUMNS | uniq | gzip -c > $TMP_NOV.tmp
    mv $TMP_NOV.tmp $TMP_NOV
fi

# Join pre and novel scored variants
{
    echo "##CADD $GENOMEBUILD-$VERSION (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health 2013-2019. All rights reserved.";
    head -n 1 < <(zcat $TMP_NOV);
    zcat $TMP_PRE $TMP_NOV | grep -v "^#" | sort -k1,1 -k2,2n -k3,3 -k4,4 || true;
} | bgzip -c > $OUTFILE;
rm $TMP_NOV
rm $TMP_PRE

OUTFILE=$(echo $OUTFILE |  sed 's/^\.\///')
echo -e "\nCADD scored variants written to file: $OUTFILE"

exit 0
