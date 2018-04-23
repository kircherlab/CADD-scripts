#!/bin/env python

import pandas as pd
from scipy.sparse import csr_matrix, vstack, save_npz
import numpy as np
import time
import sys
from argparse import ArgumentParser

def saveSparse(inputFile, outFile=None, chunksize=10000, dtype=np.float32, load=False, header=None):

    i = 0
    t0 = time.time()

    try:
        csr_block = []
        # open test data in chunks and write them to the container
        for chunk in pd.read_csv(inputFile, chunksize=chunksize, dtype=dtype, header=header):

            csr = csr_matrix(chunk)
            csr_block.append(csr)
            i += 1
            print('Read %i samples ...\n' % ((i - 1) * chunksize + chunk.shape[0]))

        t1 = time.time()
        print('Parsing took %.2f seconds\n' % (t1 - t0))

        csr = vstack(csr_block)
        del csr_block
    except pd.errors.EmptyDataError:
        sys.stderr.write('Input file is empty! Returning empty file/matrix as well.\n')
        csr = csr_matrix(np.ndarray(shape=(0,0), dtype=dtype))

    if load:
        return csr
    save_npz(outFile, csr)

if __name__ == '__main__':
    parser = ArgumentParser(description="""This script reads a csv.gz and saves it as scipy sparse matrix (*.npz)""")
    parser.add_argument("-i", "--input", dest="input", type=str, required=True,
                        help="File location of the imputed data matrix.")
    parser.add_argument("-o", "--output", dest="output", type=str,
                        required=True,
                        help="Location were the generated file is stored")
    parser.add_argument("-t", "--dtype", dest="dtype", type=str,
                        default=None,
                        help="data type the variables are stored in (float16, float64). Default is float32")
    parser.add_argument("-c", "--chunksize", dest="chunksize", type=int,
                        default=100000,
                        help="Number of lines read per chunk -> influences memory usage")
    args = parser.parse_args()

    dtype = np.float32
    if args.dtype in ['float', 'float64']:
        dtype = np.float64
    elif args.dtype == 'float16':
        dtype = np.float16
    elif args.dtype == 'float128':
        dtype = np.float128
    
    saveSparse(args.input, args.output, dtype=dtype, chunksize=args.chunksize)
