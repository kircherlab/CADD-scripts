#!/bin/env python

from argparse import ArgumentParser
from sklearn.externals import joblib
import numpy as np
import gzip
import sys
import pandas as pd

parser = ArgumentParser(description="""This script reads a list of imputed variants and predicts
                        their assignment.""")
parser.add_argument("-i", "--input", dest="input", type=str,
                    required=True, nargs='+',
                    help="File location of imputed variants as csv")
parser.add_argument("--chunksize", dest="chunksize", type=int,
                    default=100000,
                    help="Number of lines read per chunk -> influences memory usage")
parser.add_argument("-m", "--model", dest="model", type=str,
                    help="path of the model folder")
parser.add_argument("-o", "--output", dest="output", type=str,
                    default=None, nargs='+',
                    help="Location were the generated predictions are stored, default stdout")
parser.add_argument("-c", "--colname", dest="colname", default='RawScore', type=str,
                    help="the column name in the output file")
parser.add_argument("-a", "--append", dest='append', default=None, type=str,
                    nargs='+',
                    help="file to which we append the prediction at end of line")
parser.add_argument("-d", "--delimiter", dest="delimiter", default="\t",
                    type=str,
                    help="delimter of the output file if appending")
parser.add_argument("--hastarget", dest="hastarget", default=False,
                    action="store_true",
                    help="if dataset contain first column with the target variable")

args = parser.parse_args()

if args.output and len(args.input) != len(args.output):
    if len(args.output) == 1:
        sys.stderr.write('Concatenating input files\n')
    else:
        raise NotImplementedError('Cannot predict %i dataset into %i files' % (len(args.input), len(args.output)))

if args.append and len(args.input) != len(args.append):
    raise NotImplementedError('Number of input files has to be the same as the number of appending files (if the latter is >0)')

if args.output is None:
    out_file = sys.stdout
else:
    out_file = gzip.open(args.output[0], 'w')

model, scaler = joblib.load(args.model)

first_line = True

for n, infile in enumerate(args.input):

    if args.append:
        append_file = gzip.open(args.append[n], 'r')

    try:
        for chunk in pd.read_csv(infile, chunksize=args.chunksize, dtype=np.float32, header=None):

            mat_variants = np.array(chunk)

            if args.hastarget:
                mat_variants = mat_variants[:,1:]

            scaler.transform(mat_variants)

            res_variants = model.decision_function(mat_variants).flatten()

            if args.append:
                for res in res_variants:
                    while True:
                        line = append_file.readline().strip()
                        if line.startswith('#') or line.startswith('Chrom'):
                            if first_line:
                                out_file.write(''.join([line, args.delimiter,
                                                        args.colname, '\n']))
                                first_line = False
                        else:
                            break
                    out_file.write(''.join([line, args.delimiter,
                                            str(res), '\n']))

            else:
                if first_line:
                    out_file.write(args.colname + '\n')
                    first_line = False
                for res in res_variants:
                    out_file.write(str(res) + '\n')
    except pd.errors.EmptyDataError:
        sys.stderr.write('Input file %s is empty.\n' % infile)
        if first_line and args.append:
            while True:
                line = append_file.readline().strip()
                if line.startswith('#') or line.startswith('Chrom'):
                    out_file.write(''.join([line, args.delimiter,
                                            args.colname, '\n']))
                    first_line = False
                    break

    if args.output and len(args.output) > 1:
        out_file.close()
        if n + 1 < len(args.output): # open next output file if not in last iteration
            first_line = True
            out_file = gzip.open(args.output[n+1], 'w')

if args.output and len(args.output) == 1:
    out_file.close()
