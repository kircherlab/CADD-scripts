import sys, os
from argparse import ArgumentParser
from ConfigParser import ConfigParser

import vcf
import pysam
from collections import OrderedDict
import random
from lib.Annotations import annotations

# reading the situation
parser = ArgumentParser(description="%prog name")
parser.add_argument("-i", "--input", dest="input", type=str, default=None,
                    help="File with variant effect predictor output (else stdin)")
parser.add_argument("-o", "--output", dest="output", type=str, default=None,
                    help="File were to write to (else stdout)")
parser.add_argument("-c", "--config", dest="config", type=str, default=None,
                    help="File that configures the used annotation files")
parser.add_argument("--noHeader", dest="header", action="store_false",
                   default=True,
                   help="Whether to output a header of the written columns")
parser.add_argument("--continuous", dest="continuous",
                    default=False, action="store_true",
                    help="Input file is coordinate sorted and does not contain coordinate ranges -- not save for indels! (default OFF)")
parser.add_argument("--choose", dest="choose", default=False, action="store_true",
                    help="Choose a random line of the same hierarchy instead of all (def OFF)")

args = parser.parse_args()

if args.input:
    vcf_reader = vcf.Reader(open(args.input, 'r'))
else:
    vcf_reader = vcf.Reader(sys.stdin)

if args.output:
    stdout = open(args.output, 'wb')
else:
    stdout = sys.stdout

# define the used annotations:
if args.config:
    conf = ConfigParser()
    conf.read(args.config)
    pathConf = {a:b for a,b in conf.items('Path')}
    annotationConf =  {a:b for a,b in conf.items('Annotations')}

    if pathConf['roottype'].lower() in ['v', 'var', 'env']:
        root_dir = os.environ[pathConf['rootdir']]
    else:
        root_dir = pathConf['rootdir']
    root_anno = (root_dir+pathConf['annodir']).replace('//','/')
    reference = root_anno + pathConf['reference']

    included_annotations = [annotationName for annotationName in annotationConf.keys() if annotationConf[annotationName] != 'Ignore']
    for annotation in annotations:
        name = annotation.name.lower()
        if name in included_annotations and annotationConf[name] != 'True':
            annotation.path = root_anno + annotationConf[name]
else:
    # define project annotation directory
    root_dir = os.environ['CADD']
    root_anno = (root_dir+'/scoring/annotations').replace('//','/')
    for annotation in annotations:
        if hasattr(annotation, 'path'):
            annotation.path = root_anno + annotation.path

    reference = '%s/reference/reference.fa' % (root_anno)

    # use all annotations without config
    included_annotations = [annotation.name.lower() for annotation in annotations]

annotations = [annotation for annotation in annotations if annotation.mandatory or annotation.name.lower() in included_annotations]

cons_annotations = []
nocons_annotations = []
for annotation in annotations:
    if annotation.consequence:
        cons_annotations.append(annotation)
    else:
        nocons_annotations.append(annotation)

annotationFeatures = []
for annotation in annotations:
    if hasattr(annotation, 'features'):
        annotationFeatures.extend(annotation.features)
    else:
        annotationFeatures.append(annotation.name)
annotationFeatures = OrderedDict.fromkeys(annotationFeatures).keys() # remove duplicates

genome_index = pysam.FastaFile(reference)

for annotation in annotations:
    if hasattr(annotation, 'load'):
        annotation.load(args)

info_columns = vcf_reader.infos['CSQ'].desc.split('Format: ')[1].split('|')

output_columns = ['Chrom', 'Pos', 'Ref', 'Alt', 'Type']
output_columns.extend(annotationFeatures)

if args.header:
    stdout.write('#' + '\t'.join(output_columns) + '\n')

# processing
for record in vcf_reader:

    res = {}
    res['Chrom'] = record.CHROM
    if len(res['Chrom']) > 2: # quickfix: only support main chromosomes, discard others
        continue

    res['Ref'] = str(record.REF)
    res['Alt'] = str(record.ALT[0]) # there should always be only one
    res['Pos'] = record.POS
    res['Start'] = res['Pos']
    res['End'] =  res['Pos'] + len(res['Ref'])

    # sequence around the variant, so that we only retrieve it once
    # the retrieved sequence depends on the type: (uses start_seq and end_seq)
    #     SNV: pos - 1 to pos
    #     INS: pos to pos
    #     DEL: pos to pos + len(Ref)
    #     PNV: pos - 1 to pos - 1 + len(Ref)
    # the used segment for tabix calls is slightly different: (uses res['Start'] and res['End'])
    #     SNV: pos
    #     PNV: pos to pos + len(Ref) - 1
    #     INS: pos to pos + 1 (the two bases around the insertion)
    #     DEL: pos + 1 to pos + len(Ref) - 1
    # complex INS are like DEL
    if len(res['Ref']) == len(res['Alt']):
        res['Type'] = 'SNV'
        start_seq = res['Pos'] - 1
        end_seq = res['Pos'] + len(res['Ref']) - 1
        res['Start'] = res['Pos']
        res['End'] = end_seq
    else:
        start_seq = res['Pos']
        if len(res['Ref']) == 1:
            res['Type'] = 'INS'
            end_seq = res['Pos']
            res['Start'] = res['Pos']
            res['End'] = res['Pos'] + 1
        else:
            if len(res['Alt']) == 1 or len(res['Ref']) > len(res['Alt']):
                res['Type'] = 'DEL'
            else: # Exceptions are complex INS like AG -> CTGCT (same with DEL in previous if)
                res['Type'] = 'INS'
            end_seq = res['Pos'] + len(res['Ref']) - 1
            res['Start'] = res['Pos'] + 1
            res['End'] = end_seq
    try:
        if start_seq > 74:
            res['Seq'] = genome_index.fetch(res['Chrom'], start_seq-75, end_seq+75).upper()
            if len(res['Seq']) < 75:
                sys.stderr.write('Encountered variant outside of chromosome boundary: %s\t%i\t%s\t%s\n' % (res['Chrom'], res['Pos'], res['Ref'], res['Alt']))
                continue
        else:
            res['Seq'] = (75-start_seq) * 'N' + genome_index.fetch(res['Chrom'], 0, end_seq+75).upper()
    except KeyError:
        sys.stderr.write('Encountered unknown chromosome name %s: %s\t%i\t%s\t%s\n' % (res['Chrom'],res['Chrom'], res['Pos'], res['Ref'], res['Alt']))
        continue

    # enable annotations to access the vcf info column
    res['Info'] = record.INFO

    # read those annotations, that are independent of the consequence
    for annotation in nocons_annotations:
        res = annotation.process(res)

    # evalute the different annotated consequences separately
    org_res = res
    res_list = []
    for csq in record.INFO['CSQ']:
        res = org_res.copy()

        # read info fields
        info = csq.split("|")
        for num, col in enumerate(info_columns):
            res[col] = info[num]

        # read the consequence specific annotation annotations
        for annotation in cons_annotations:
            res = annotation.process(res)

        res_list.append(res)

    # add motif features from separate consequence lines
    num_motiffeature = 0
    motif_res = None
    for res in res_list:
        if res['Feature_type'] == 'MotifFeature':
            num_motiffeature += 1
            if motif_res and 'motifEScoreChng' in motif_res: # proliferate the highest scoring motif
                if 'motifEScoreChng' in res and res['motifEScoreChng'] > motif_res['motifEScoreChng']:
                    motif_res = res
            else:
                motif_res = res
    if 0 < num_motiffeature < len(res_list): # ignore/remove motif features if not all consequences are such
        res_list2 = []
        transfers = [t for t in ['motifEName','motifEHIPos','motifEScoreChng'] if t in motif_res]
        for res in res_list:
            if res['Feature_type'] != 'MotifFeature':
                for t in transfers:
                    res[t] = motif_res[t]
                    res['motifECount'] = num_motiffeature
                res_list2.append(res)
        res_list = res_list2

    # when choosing only one functional annotation per variant
    if args.choose:
        consequence_levels = [['STOP_GAINED','STOP_LOST','FRAME_SHIFT',
                               'INFRAME', 'NON_SYNONYMOUS',
                               'NONCODING_CHANGE', 'CANONICAL_SPLICE',
                               'SPLICE_SITE', 'SYNONYMOUS', 'UNKNOWN'],
                              ['5PRIME_UTR', '3PRIME_UTR', 'REGULATORY'],
                              ['INTRONIC', 'UPSTREAM', 'DOWNSTREAM']]
        res_list2 = []
        top_level = 3
        for res in res_list:
            for level, conseqs in enumerate(consequence_levels):
                if level > top_level:
                    break
                if res['Consequence'] in conseqs:
                    if level == top_level:
                        res_list2.append(res)
                    else:
                        top_level = level
                        res_list2 = [res]
                    break
        if top_level < 3:
            res_list = res_list2
        res_list = random.sample(res_list, 1)

    # write off to stream
    for res in res_list:
        res = {r: str(val) for r, val in res.items()}
        res = {r: val.rstrip('0').rstrip('.') if '.0' in val else val for r, val in res.items()}
        stdout.write('\t'.join([res[col] if col in res else "NA" for col in output_columns]) + '\n')

# cleaning up
genome_index.close()
if args.output:
    stdout.close()
