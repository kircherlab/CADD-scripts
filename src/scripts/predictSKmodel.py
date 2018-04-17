#!/bin/env python

from argparse import ArgumentParser
from sklearn.externals import joblib
from scipy.sparse import load_npz
import gzip
import sys

parser = ArgumentParser(description="""This script reads a list of imputed variants and predicts
                        their assignment.""")
parser.add_argument("-i", "--input", dest="input", type=str,
                    required=True, nargs='+',
                    help="File location of imputed variants as npz (sparse or dense)")
parser.add_argument("-m", "--model", dest="model", type=str,
                    help="path of the model folder")
parser.add_argument("-o", "--output", dest="output", type=str,
                    default=None, nargs='+',
                    help="Location were the generated predictions are stored, default stdout")
parser.add_argument("-a", "--append", dest='append', default=None, type=str,
                    nargs='+',
                    help="file to which we append the prediction at end of line")
parser.add_argument("-d", "--delimiter", dest="delimiter", default="\t",
                    type=str,
                    help="delimter of the output file if appending")

args = parser.parse_args()

if args.output and len(args.input) != len(args.output):
    if len(args.output) == 1:
        print('Concatenating input files')
    else:
        raise NotImplementedError('Cannot predict %i dataset into %i files' % (len(args.input), len(args.output)))

if args.append and len(args.input) != len(args.append):
    raise NotImplementedError('Number of input files has to be the same as the number of appending files (if the latter us >0)')

if args.output is None:
    out_file = sys.stdout
else:
    out_file = gzip.open(args.output[0], 'w')

for n, infile in enumerate(args.input):

    mat_variants = load_npz(infile).toarray()

    model, scaler = joblib.load(args.model)

    scaler.transform(mat_variants)

    try:
        res_variants = model.predict_proba(mat_variants)[:,1]
    except AttributeError:
        res_variants = model.decision_function(mat_variants).flatten()

    num = 0
    if args.append:
        with gzip.open(args.append[n], 'r') as in_file:
            for line in in_file:
                line = line.strip()
                if line.startswith('#') or line.startswith('Chrom'):
                    if n == 0:
                        out_file.write(''.join([line, args.delimiter,
                                                'RawScore', '\n']))
                else:
                    out_file.write(''.join([line, args.delimiter,
                                            str(res_variants[num]), '\n']))
                    num += 1
    else:
        out_file.write('RawScore' + '\n')
        for i in res_variants:
            out_file.write(str(i) + '\n')

    if args.output and len(args.output) > 1:
        out_file.close()
        if n + 1 < len(args.output): # open next output file if not in last iteration
            out_file = gzip.open(args.output[n+1], 'w')
