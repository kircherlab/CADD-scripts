#!/bin/bash

usage="$(basename "$0") [-o <outfile>] [-g <genomebuild>] [-v <caddversion>] [-a] <infile>  -- CADD version 1.6

where:
    -h  show this help text
    -o  out tsv.gz file (generated from input file name if not set)
    -g  genome build (supported are GRCh37 and GRCh38 [default: GRCh38])
    -v  CADD version (only v1.6 possible with this set of scripts [default: v1.6])
    -a  include annotation in output
        input vcf of vcf.gz file (required)
    -q  print basic information about snakemake run
    -p  print full information about the snakemake run
    -c  number of cores that snakemake is allowed to use [default: 1]
    "

unset OPTARG
unset OPTIND
export LC_ALL=C

GENOMEBUILD="GRCh38"
ANNOTATION=false
OUTFILE=""
VERSION="v1.6"
VERBOSE="-q"
CORES="1"
while getopts ':ho:g:v:c:aqp' option; do
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
    v) CORES=$OPTARG
       ;;
    a) ANNOTATION=true
       ;;
    q) VERBOSE=""
       ;;
    p) VERBOSE="-p"
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND-1))

INFILE=$1

echo "CADD-v1.6 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health 2013-2020. All rights reserved."

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

if [ "$VERSION" != "v1.6" ]
then
    echo "Unknown/Unsupported CADD version $VERSION. This set of script currently only supports v1.6."
    echo "If you want to score another version of CADD, please download the accordingly tagged version of the scripts"
    exit 1
fi

if [ "$ANNOTATION" = 'true' ]
then
    CONFIG=$CADD/config/config_${GENOMEBUILD}_${VERSION}.yml
else
    CONFIG=$CADD/config/config_${GENOMEBUILD}_${VERSION}_noanno.yml
fi

# Setup temporary folder that is removed reliably on exit and is outside of
# the CADD-scripts directory.
TMP_FOLDER=$(mktemp -d)
trap "rm -rf $TMP_FOLDER" ERR EXIT

# Temp files
TMP_INFILE=$TMP_FOLDER/$NAME.$FILEFORMAT
TMP_OUTFILE=$TMP_FOLDER/$NAME.tsv.gz

cp $INFILE $TMP_INFILE

echo "Running snakemake pipeline:"
echo snakemake $TMP_OUTFILE --use-conda --conda-prefix $CADD/envs --cores $CORES
echo --configfile $CONFIG --snakefile $CADD/Snakefile $VERBOSE
snakemake $TMP_OUTFILE --use-conda --conda-prefix $CADD/envs --cores $CORES \
    --configfile $CONFIG --snakefile $CADD/Snakefile $VERBOSE

mv $TMP_OUTFILE $OUTFILE
rm $TMP_INFILE # is in temp folder, should not be necessary

OUTFILE=$(echo $OUTFILE |  sed 's/^\.\///')
echo -e "\nCADD scored variants written to file: $OUTFILE"

exit 0
