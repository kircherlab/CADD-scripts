#!/usr/bin/env python

from ConfigParser import ConfigParser
from argparse import ArgumentParser
from lib.columnInfo import get_columns
import sys
import numpy as np
import random
import gzip

parser = ArgumentParser(description="%prog name")
parser.add_argument("-i", "--input", dest="input", type=str,
                    default=None,
                    help="File location of variant data in tsv format. If not specified use stdin")
parser.add_argument("-o", "--output", dest="output", type=str,
                    default=None,
                    help="Location were the generated file is stored. If not specified use stdout")
parser.add_argument("-c", "--config", dest="config", type=str,
                    default='lib/exampleTracks.cfg',
                    help="Config file that specifies used tracks")
parser.add_argument("-b", "--cat2bool", dest="cat2bool", action='store_true',
                    help="Specify whether categories are split into multiple boolean classifier")
parser.add_argument("--noheader", dest="noheader", default=False,
                    action="store_true",
                    help="Do not print header line")

args = parser.parse_args()

# define input and output sources
if args.input is None:
    stdin = sys.stdin
else:
    stdin = gzip.open(args.input, 'r')

if args.output is None:
    stdout = sys.stdout
else:
    stdout = gzip.open(args.output, 'w')

### TRACK PREPARATION ###
newColumnNames, _, trackData, config = \
    get_columns(args.config, None, args.cat2bool, False)

### DATA PROCESSING ###
columnNames = []
for line in stdin:
    line = line.strip()
    # detect header line
    if line.startswith('#') or line.startswith('Chrom'):
        columnNames = line.strip('#').split('\t')
        columnNames = np.array(map(str.lower, columnNames))
        if not args.noheader:
            stdout.write(','.join(newColumnNames) + '\n')
        continue
    # associate fields with column names
    fieldsDict = dict(zip(columnNames, line.split('\t')))
    outFields = []
    indicatorFields = []

    for trackName, status in config['Tracks']:
        if 'colname' in trackData[trackName].keys():
            trackName = trackData[trackName]['colname']
        track = trackData[trackName]

        if status == 'Ignore':
            # outFields.append(fieldsDict[trackName])
            continue
        if status != 'True':
            continue

        if track['type'] == 'combined':
            if args.cat2bool:
                i = trackData[track['base']]['id']
                baseArray = np.array(outFields[i:i+len(trackData[track['base']]['categories'])])
            else:
                baseValue = outFields[trackData[track['base']]['id']]
                baseArray = (np.array(trackData[track['base']]['categories']) == baseValue).astype(int)
            values = []
            for child in track['child']:
                values.extend(baseArray * outFields[trackData[child]['id']])
            outFields.extend([value for value in values])
        else:
            try:
                if 'derive' in track.keys():
                    value = track['derive'](fieldsDict)
                else:
                    value = fieldsDict[trackName]
                if track['type'] in [float, int]:
                    value = track['type'](value)
                if 'transformation' in track.keys():  # transform is slightly redundant to derive
                    value = track['transformation'](value)
                if track['type'] is list:
                    assert (value in track['categories'])
                if 'indicator' in track.keys():
                    indicatorFields.append('0')
            except:
                value = track['na_value']
                if 'indicator' in track.keys():
                    indicatorFields.append('1')

            if args.cat2bool and track['type'] is list:
                values = (np.array(track['categories']) == value).astype(int)
                outFields.extend(values)
            else:
                outFields.append(value)

    # minimize zeros and stringify
    outFields = ['0' if f == 0 else str(f) for f in outFields]

    outFields.extend(indicatorFields)

    stdout.write(','.join(outFields) + '\n')

if args.input is not None:
    stdin.close()

elif args.output is not None:
    stdout.close()
