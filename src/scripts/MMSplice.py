#!/bin/env python3

# Import
from tqdm import tqdm
from concise.preprocessing import encodeDNA
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice
from mmsplice.utils import logit, predict_deltaLogitPsi, \
    predict_pathogenicity, predict_splicing_efficiency
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

def predict_batch_fast(model, dataloader, batch_size=512, progress=True,
                       splicing_efficiency=False):
    """
    Return the prediction as a table
    Args:
      model: mmsplice model object.
      dataloader: dataloader object.
      progress: show progress bar.
      splicing_efficiency: adds splicing_efficiency prediction as column
    Returns:
      iterator of pd.DataFrame of modular prediction, delta_logit_psi,
        splicing_efficiency, pathogenicity.
    """
    dataloader.encode = False
    dt_iter = dataloader.batch_iter(batch_size=batch_size)
    if progress:
        dt_iter = tqdm(dt_iter)

    ref_cols = ['ref_acceptorIntron', 'ref_acceptor',
                'ref_exon', 'ref_donor', 'ref_donorIntron']
    alt_cols = ['alt_acceptorIntron', 'alt_acceptor',
                'alt_exon', 'alt_donor', 'alt_donorIntron']

    cat_list = ['acceptor_intron', 'acceptor', 'exon', 'donor', 'donor_intron']
    cats = {'acceptor_intron': lambda x: model.acceptor_intronM.predict(x),
            'acceptor': lambda x: logit(model.acceptorM.predict(x)),
            'exon': lambda x: model.exonM.predict(x),
            'donor': lambda x: logit(model.donorM.predict(x)),
            'donor_intron': lambda x: model.donor_intronM.predict(x)}

    for batch in dt_iter:
        refs,alts = {},{}
        for cat, model_eval in cats.items():
            alterations = batch['inputs']['seq'][cat] != \
                          batch['inputs']['mut_seq'][cat]
            if np.any(alterations):
                sequences = list(set([str(s) for s in list(batch['inputs']['seq'][cat][alterations]) + \
                      list(batch['inputs']['mut_seq'][cat][alterations])]))
                prediction = model_eval(encodeDNA(sequences)).flatten()
                pred_dict = {s: p for s, p in zip(sequences, prediction)}

                refs[cat] = [pred_dict[batch['inputs']['seq'][cat][i]] \
                            if a else 0 for i, a in enumerate(alterations)]
                alts[cat] = [pred_dict[batch['inputs']['mut_seq'][cat][i]] \
                            if a else 0 for i, a in enumerate(alterations)]
            else:
                refs[cat] = [0] * len(alterations)
                alts[cat] = [0] * len(alterations)
        X_ref = np.array([refs[cat] for cat in cat_list]).T
        X_alt = np.array([alts[cat] for cat in cat_list]).T
        ref_pred = pd.DataFrame(X_ref, columns=ref_cols)
        alt_pred = pd.DataFrame(X_alt, columns=alt_cols)

        df = pd.DataFrame({
            'ID': batch['metadata']['variant']['STR'],
            'exons': batch['metadata']['exon']['annotation'],
        })
        for k in ['exon_id', 'gene_id', 'gene_name', 'transcript_id']:
            if k in batch['metadata']['exon']:
                df[k] = batch['metadata']['exon'][k]

        df['delta_logit_psi'] = predict_deltaLogitPsi(X_ref, X_alt)
        df = pd.concat([df, ref_pred, alt_pred], axis=1)

        # pathogenicity does not work
        #if pathogenicity:
        #    df['pathogenicity'] = predict_pathogenicity(X_ref, X_alt)

        if splicing_efficiency:
            df['efficiency'] = predict_splicing_efficiency(X_ref, X_alt)

        yield df

def predict_table_fast(model,
                      dataloader,
                      batch_size=512,
                      progress=True,
                      pathogenicity=False,
                      splicing_efficiency=False):
    """
    Return the prediction as a table
    Args:
      model: mmsplice model object.
      dataloader: dataloader object.
      progress: show progress bar.
      splicing_efficiency: adds  splicing_efficiency prediction as column
    Returns:
      pd.DataFrame of modular prediction, delta_logit_psi, splicing_efficiency,
        pathogenicity.
    """
    return pd.concat(predict_batch_fast(model, dataloader, batch_size=batch_size,
                                   progress=progress,
                                   splicing_efficiency=splicing_efficiency))

parser = ArgumentParser(description="%prog name")
parser.add_argument("-o", "--output", dest="output", type=str,
                    help="Output vcf file {default stdout}")
parser.add_argument("-i", "--input", dest="input", type=str, required=True,
                    help="Input vcf.gz file")
parser.add_argument("-f", "--fasta", dest="fasta", type=str, required=True,
                    help="Genome fasta")
parser.add_argument("-g" "--gtf", dest="gtf", type=str, required=True,
                    help="GTF file with all exons (may be the pickled version output by the script)")
parser.add_argument("-d", "--detail", dest="detail", action='store_true',
                    help="Do not report submodel score differences but delta_logit_psi and ref and alt scores separately")

args = parser.parse_args()

dl = SplicingVCFDataloader(args.gtf,
                           args.fasta,
                           args.input)

# Specify model
model = MMSplice()

try:
    # Do prediction
    predictions = predict_table_fast(model, dl, batch_size=512)

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
    if args.detail:
        for line in reader:
            if line.startswith("#"):
                if line.startswith('#C'):
                    writer.write('##INFO=<ID=MMSplice,Number=.,Type=String,Description="MMSplice scores. These include scores for ref (R) and alt (A) separately. Format: SYMBOL|delta_logit_psi|R_acceptorIntron|R_acceptor|R_exon|R_donor|R_donorIntron|A_acceptorIntron|A_acceptor|A_exon|A_donor|A_donorIntron">\n')

                writer.write(line)
                continue
            fields = line.strip().split('\t')
            var = "%s:%s:%s:['%s']" % (fields[0], fields[1], fields[3], fields[4])
            if var in pred_dict:
                val_list = []
                for p in pred_dict[var]:
                    val_list.append('|'.join([p[0]] + ['%.3f' % p for p in p[1:]]))
                val = ','.join(val_list)
                if len(fields) >= 8:
                    fields[7] = fields[7] + ';MMSplice=' + val
                else:
                    fields.extend(['.']*(7 - len(fields)))
                    fields.append('MMSplice=' + val)
            writer.write('\t'.join(fields) + '\n')
    else:
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
