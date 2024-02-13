#!/bin/env python3

# Import
from tqdm import tqdm
# from concise.preprocessing import encodeDNA
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table
import pandas as pd
import numpy as np

from argparse import ArgumentParser
import sys
import gzip

def max_geneEff(df):
    """ Summarize largest absolute effect per variant per gene across all affected exons.
    Similar to mmsplice.utils.max_varEff
    Args:
        df: result of `predict_all_table`
    """
    df_max = df.groupby(['ID', 'gene_name'], as_index=False).agg(
        {'delta_logit_psi': lambda x: max(x, key=abs)})

    df_max = df_max.merge(df, how='left', on=['ID', 'gene_name', 'delta_logit_psi'])
    df_max = df_max.drop_duplicates(subset=['ID', 'gene_name', 'delta_logit_psi'])
    return df_max

parser = ArgumentParser(description="%prog name")
parser.add_argument("-o", "--output", dest="output", type=str,
                    help="Output vcf file {default stdout}")
parser.add_argument("-i", "--input", dest="input", type=str, required=True,
                    help="Input vcf.gz file")
parser.add_argument("-f", "--fasta", dest="fasta", type=str, required=True,
                    help="Genome fasta")
parser.add_argument("-g" "--gtf", dest="gtf", type=str, required=True,
                    help="GTF file with all exons (may be the pickled version output by the script)")

args = parser.parse_args()

dl = SplicingVCFDataloader(args.gtf,
                           args.fasta,
                           args.input)

# Specify model
model = MMSplice()

try:
    # Do prediction
    predictions = predict_all_table(model, dl, batch_size=512)

    # Summerize with maximum effect size
    predictionsMax = max_geneEff(predictions)

    pred_dict = {}
    for p in predictionsMax[['ID', 'gene_name', 'delta_logit_psi', 'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron']].values:
        if p[0] in pred_dict:
            pred_dict[p[0]].append(p[1:])
        else:
            pred_dict[p[0]] = [p[1:]]
except ValueError: # it can happen that no variant is in a splice region an therefore recieves a score leading to an empty array
    pred_dict = {}

if args.output:
    writer = gzip.open(args.output, 'wt')
else:
    writer = sys.stdout

# write prediction to stdout or file
with gzip.open(args.input, 'rt') as reader:
    for line in reader:
        if line.startswith("#"):
            if line.startswith('#C'):
                writer.write('##INFO=<ID=MMSplice,Number=.,Type=String,Description="MMSplice scores for five submodels. Format: SYMBOL|acceptorIntron|acceptor|exon|donor|donorIntron">\n')

            writer.write(line)
            continue
        fields = line.strip().split('\t')
        var = "%s:%s:%s:['%s']" % (fields[0], fields[1], fields[3], fields[4])
        if var in pred_dict:
            val_list = []
            for p in pred_dict[var]:
                scores = [p[i+7] - p[i+2] for i in range(5)]
                val_list.append('|'.join([p[0]] + ['%.3f' % s for s in scores]))
            val = ','.join(val_list)
            if len(fields) >= 8:
                fields[7] = fields[7] + ';MMSplice=' + val
            else:
                fields.extend(['.']*(7 - len(fields)))
                fields.append('MMSplice=' + val)
        writer.write('\t'.join(fields) + '\n')
