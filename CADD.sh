#!/bin/bash

usage="$(basename "$0") [-o <outfile>] [-g <genomebuild>] [-v <caddversion>] [-a] <infile>  -- CADD version 1.7

where:
    -h  show this help text
    -o  out tsv.gz file (generated from input file name if not set)
    -g  genome build (supported are GRCh37 and GRCh38 [default: GRCh38])
    -v  CADD version (only v1.7 possible with this set of scripts [default: v1.7])
    -a  include annotation in output
        input vcf of vcf.gz file (required)
    -m  use conda only (no apptainer/singularity)
    -r  singularity/apptainer arguments, e.g. \"--bind /data/mnt/x --nv\" [default \"\" but will always add \"--bind $TMP_DIR\"]
    -q  print basic information about snakemake run
    -p  print full information about the snakemake run
    -d  do not remove temporary directory for debug puroposes
    -c  number of cores that snakemake is allowed to use [default: 1]
    "

unset OPTARG
unset OPTIND
export LC_ALL=C

GENOMEBUILD="GRCh38"
ANNOTATION=false
CONDAONLY=false
OUTFILE=""
VERSION="v1.7"
SIGNULARITYARGS=""
VERBOSE="-q"
CORES="1"
RM_TMP_DIR=true
while getopts ':ho:g:v:c:amr:qpd' option; do
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
    c) CORES=$OPTARG
       ;;
    a) ANNOTATION=true
       ;;
    m) CONDAONLY=true
       ;;
    r) SIGNULARITYARGS=$OPTARG
       ;;
    q) VERBOSE=""
       ;;
    p) VERBOSE="-p"
       ;;
    d) RM_TMP_DIR=false
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND-1))

INFILE=$1

echo "CADD-v1.7 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health at Charite - Universitatsmedizin Berlin 2013-2024. All rights reserved."

set -ueo pipefail

# check if input file does exist
if [ "$INFILE" == "" ]
then
    echo "No input file specified. To run CADD, a list of variants has to be provided in a vcf or vcf.gz file."
    exit 1
elif [ ! -f "$INFILE" ]
then
    echo "Input file $INFILE does not exist."
    exit 1
fi

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

if [ "$VERSION" != "v1.7" ]
then
    echo "Unknown/Unsupported CADD version $VERSION. This set of script currently only supports v1.7."
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
if [ "$RM_TMP_DIR" = 'true' ]
then
    trap "rm -rf $TMP_FOLDER" ERR EXIT
fi

# Temp files
TMP_INFILE=$TMP_FOLDER/$NAME.$FILEFORMAT
TMP_OUTFILE=$TMP_FOLDER/$NAME.tsv.gz

cp $INFILE $TMP_INFILE

# setup bindings of singularity args
if [ "$CONDAONLY" = 'true' ]
then
    SIGNULARITYARGS=""
else
    SIGNULARITYARGS="apptainer --apptainer-prefix $CADD/envs/apptainer --singularity-args \"--bind ${TMP_FOLDER} ${SIGNULARITYARGS}\""
fi

echo "Running snakemake pipeline:"

command="snakemake $TMP_OUTFILE \
    --sdm conda $SIGNULARITYARGS --conda-prefix $CADD/envs/conda \
    --cores $CORES --configfile $CONFIG \
    --snakefile $CADD/Snakefile $VERBOSE"

echo -e $command

eval $command

mv $TMP_OUTFILE $OUTFILE
rm $TMP_INFILE # is in temp folder, should not be necessary

OUTFILE=$(echo $OUTFILE |  sed 's/^\.\///')
echo -e "\nCADD scored variants written to file: $OUTFILE"

exit 0
