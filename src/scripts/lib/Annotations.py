import os, sys
import pysam
import numpy as np
from AnalysisLib import get_from_tabix,get_range_from_tabix
from AnalysisLib import eval_dbSNP_frequencies,eval_ESP_frequencies

###################
### ANNOTATIONS ###
###################

class Annotation(object):
    consequence = False
    mandatory = False

class FeatureAnnotation(Annotation):

    def process(self, res):
        self._retrieve(res)
        return self._get_score(res)

class TabixAnnotation(FeatureAnnotation):
    multirange = False
    rangescore = False
    zerobased = False

    def load(self, args):
        if self.path != '' and (not os.path.exists(self.path) or not os.path.exists(self.path+".tbi")):
            sys.stderr.write("%s annotation: Require valid path to compressed tabix file and tabix index file.\n" % self.name)
        elif os.path.exists(self.path) and os.path.exists(self.path+".tbi"):
            self.tabix = pysam.TabixFile(self.path,'r'),None,None,None,None,None,None,"%s annotation" % self.name
            self.continuous = args.continuous

    def _retrieve(self, res):
        self.tabix = get_range_from_tabix(self.tabix,
                                          res['Chrom'],
                                          res['Start'],
                                          res['End'],
                                          rangescore=self.rangescore,
                                          continuous=self.continuous,
                                          multirange=self.multirange,
                                          zerobased=self.zerobased)
        self.score = self.tabix[6]

class PosTabixAnnotation(TabixAnnotation):
    # Annotation for this exact position and not the variants whole range

    def _retrieve(self, res):
        self.tabix = get_from_tabix(self.tabix,
                                    res['Chrom'],
                                    res['Pos'],
                                    continuous=self.continuous,
                                    zerobased=self.zerobased)
        self.score = self.tabix[6]

def maxL(valuelist):
    if len(valuelist) > 0: return max(valuelist)
    else: return None

def minL(valuelist):
    if len(valuelist) > 0: return min(valuelist)
    else: return None



class EncodeAnnotation(FeatureAnnotation):

    def load(self, args):
        self.tabixlist = []
        if os.path.exists(self.path):
            cpath = "/".join(self.path.split('/')[:-1])
            infile = open(self.path)
            for filename in infile:
                fullname = cpath+'/'+filename.strip()
                if os.path.exists(fullname):
                    self.tabixlist.append((pysam.TabixFile(fullname,'r'),None,None,None,None,None,None,"%s %s"%(self.name, fullname)))
                else:
                    sys.stderr.write("Not found: %s\n"%fullname)
            infile.close()
        else:
            sys.stderr.write("Index file not found: %s\n"%self.path)
        self.continuous = args.continuous

    def _retrieve(self, res):
        newlist = []
        values = []
        for encodeTabix in self.tabixlist:
            encodeTabix = get_range_from_tabix(encodeTabix,
                                               res['Chrom'],
                                               res['Start'],
                                               res['End'],
                                               rangescore=True,
                                               continuous=self.continuous)
            encode_out = encodeTabix[6]
            if encode_out:
                values.append(maxL(map(lambda x:float(x[-1]),encode_out)))
            newlist.append(encodeTabix)
        self.tabixlist, self.score = newlist, maxL(values)

    def _get_score(self, res):
        if not self.score is None:
            res[self.name] = self.score # "%.2f" % self.score in old version
        return res

class ScoreLowest():
    def _get_score(self, res):
        if len(self.score) == 1: res[self.name] = self.score[0][-1]
        elif len(self.score) > 1:
            helper = map(lambda x:(float(x[-1]),x[-1]),self.score)
            helper.sort()
            res[self.name] = helper[0][-1]
        return res

class ScoreHighest():
    def _get_score(self, res):
        if len(self.score) == 1: res[self.name] = self.score[0][-1]
        elif len(self.score) > 1:
            helper = map(lambda x:(float(x[-1]),x[-1]),self.score)
            helper.sort()
            res[self.name] = helper[-1][-1]
        return res

class ScoreHighestPerFeature():
    datatype = float
    def _get_score(self, res):
        for i, feature in enumerate(reversed(self.features), start=1):
            if len(self.score) == 1:
                res[feature] = self.score[0][-i]
            elif len(self.score) > 1:
                scores = [self.datatype(sc[-i]) for sc in self.score if sc[-i] != 'NA']
                if scores:
                    res[feature] = str(max(scores))
        return res

class ScoreAveragePerFeature():
    def _get_score(self, res):
        for i, feature in enumerate(reversed(self.features), start=1):
            if len(self.score) == 1:
                res[feature] = self.score[0][-i]
            elif len(self.score) > 1:
                scores = [float(sc[-i]) for sc in self.score if sc[-i] != 'NA']
                if scores:
                    res[feature] = str(np.mean(scores))
        return res

class ScoreHit():
    def _get_score(self, res):
        if self.score:
            res[self.name] = 1
        return res

class AACombination(FeatureAnnotation):

    def load(self, args):
        self.combinations = {}
        if os.path.exists(self.path):
            infile = open(self.path)
            for line in infile:
                if line.startswith('#'):
                    continue
                fields = line.split()
                if len(fields) == 2:
                    AAs = tuple(fields[0].upper().split('-'))
                    self.combinations[AAs]=fields[1]
                else:
                    sys.stderr.write('%s scores, unexpected line, skipping: %s\n'%(self.name, line.strip()))
            infile.close()
        else:
            sys.stderr.write('%s scores input file does not exist: %s\n'%(self.name, self.path))

    def process(self, res):
        if 'oAA' in res and (res['oAA'],res['nAA']) in self.combinations:
            res[self.name] = str(self.combinations[(res['oAA'],res['nAA'])])
        return res

###########################
### FEATURE ANNOTATIONS ###
###########################

### REDUNDANT ANNOTATION, later derived in tracks
transitions = [('C','T'),('T','C'),('G','A'),('A','G')]
class Transversion(Annotation):
    name = 'isTv'

    def process(self, res):
        if res['Type'] == 'SNV' and len(res['Ref']) == 1:
            if (res['Ref'], res['Alt']) in transitions:
                res['isTv'] = 'FALSE'
            else:
                res['isTv'] = 'TRUE'
        return res

class Length(Annotation):
    name = 'Length'

    def process(self, res):
        res['Length']=str(max(len(res['Ref']),len(res['Alt']))-1)
        return res

noncoding_transcript_types = set(['mature_miRNA_variant',
                                 'non_coding_transcript_exon_variant'])
coding_transcript_types = set(['missense_variant',
                               'synonymous_variant',
                               'stop_gained',
                               'stop_lost',
                               'initiator_codon_variant',
                               'stop_retained_variant',
                               'frameshift_variant',
                               'inframe_insertion',
                               'inframe_deletion',
                               'incomplete_terminal_codon_variant',
                               'protein_altering_variant'
                               ])
consequence_types = [ # category, [consequences], ConsScore
                     ('STOP_GAINED', ['stop_gained'], '8'),
                     ('NON_SYNONYMOUS', ['missense_variant', 'start_lost'], '7'),
                     ('FRAME_SHIFT', ['frameshift_variant','protein_altering_variant'], '7'),
                     ('STOP_LOST', ['stop_lost', 'incomplete_terminal_codon_variant'], '7'),
                     ('INFRAME', ['inframe_insertion', 'inframe_deletion'], '6'),
                     ('CANONICAL_SPLICE', ['splice_acceptor_variant', 'splice_donor_variant'], '6'),
                     ('SPLICE_SITE', ['splice_region_variant'], '5'),
                     ('NONCODING_CHANGE', ['non_coding_transcript_exon_variant', 'mature_miRNA_variant'], '5'),
                     ('SYNONYMOUS', ['synonymous_variant', 'stop_retained_variant'], '5'),
                     ('REGULATORY', ['regulatory_region_variant', 'TF_binding_site_variant'], '4'),
                     ('5PRIME_UTR', ['5_prime_UTR_variant'], '3'),
                     ('3PRIME_UTR', ['3_prime_UTR_variant'], '2'),
                     ('INTRONIC', ['intron_variant'], '2'),
                     ('NONCODING_CHANGE', ['non_coding_transcript_variant'], '5'),
                     ('UPSTREAM', ['upstream_gene_variant'], '1'),
                     ('DOWNSTREAM', ['downstream_gene_variant'], '1'),
                     ('UNKNOWN', ['coding_sequence_variant'], '5'),
                     ('INTERGENIC', ['intergenic_variant'], '0'),
                     ]
class Consequence(Annotation):
    name = 'Consequence'
    features = ['AnnoType', 'Consequence', 'ConsScore', 'ConsDetail']
    consequence = True

    def process(self, res):
        conseqs = set(res['Consequence'].split('&'))

        details = []
        for detail in conseqs:
            for rep in ['_variant', '_region', '_gene', '_transcript']:
                detail = detail.replace(rep, '')
            details.append(detail)
        res['ConsDetail'] = ','.join(details)

        res['AnnoType'] = res['Feature_type']
        if res['AnnoType'] == 'Transcript':
            if not conseqs.isdisjoint(noncoding_transcript_types):
                res['AnnoType'] = 'NonCodingTranscript'
            elif not conseqs.isdisjoint(coding_transcript_types):
                res['AnnoType'] = 'CodingTranscript'

        for cat, cons, score in consequence_types:
            if not conseqs.isdisjoint(cons):
                res['Consequence'] = cat
                res['ConsScore'] = score
                break

        if res['AnnoType'] == '':
            res['AnnoType'] = 'Intergenic'
            res['ConsScore'] = '0'
        if res['Consequence'] in ['UPSTREAM','DOWNSTREAM']:
            res['AnnoType'] = 'Intergenic'

        return res

class GC_CpG(Annotation):
    name = 'GC_CpG'
    features = ['GC', 'CpG']

    def process(self, res):
        seq_len = float(len(res['Seq']))
        num_n = res['Seq'].count('N')
        res['GC'] = (res['Seq'].count('C') + res['Seq'].count('G') + num_n*0.41) / seq_len
        res['CpG'] = (res['Seq'].count('CG') + num_n * 0.01) / (seq_len - 1) * 2

        # reduce precision
        res['GC'] = ('%.3f' % res['GC']).rstrip('0').rstrip('.')
        res['CpG'] = ('%.3f' % res['CpG']).rstrip('0').rstrip('.')
        return res

class oldGC_CpG(GC_CpG):
    name = 'oldGC_CpG'

    def process(self, res):
        seq_len = float(len(res['Seq']))
        num_n = res['Seq'].count('N')
        res['GC'] = (res['Seq'].count('C') + res['Seq'].count('G') + num_n*0.41) / seq_len
        res['CpG'] = (res['Seq'].count('CG') + num_n * 0.01) / (seq_len - 1) * 2
        return res

class MotifScores(Annotation):
    name = 'MotifScores'
    features = ['motifECount','motifEName','motifEHIPos','motifEScoreChng']
    consequence = True

    def process(self, res):
        if res['MOTIF_NAME'] != '':
            res['motifECount'] = '1' # alternative for companyless motif features
            res['motifEName'] = res['MOTIF_NAME']
        if res['HIGH_INF_POS'] != '':
            res['motifEHIPos'] = 1 if "Y" == res['HIGH_INF_POS'] else 0
        if res['MOTIF_SCORE_CHANGE'] != '':
            res['motifEScoreChng'] = res['MOTIF_SCORE_CHANGE']
        return res

class AminoAcids(Annotation):
    name = 'AminoAcids'
    features = ['oAA', 'nAA']
    consequence = True
    mandatory=True

    def process(self, res):
        if res['Amino_acids'] != '':
            aa = res['Amino_acids'].split('/')
            if len(aa) == 2:
                res['oAA'], res['nAA'] = aa
            else:
                res['oAA'], res['nAA'] = aa[0], aa[0]
        return res

class VEPAnnotations(Annotation):
    '''This is mostly a placeholde annotation that adds annotations from the VEP
    output as features for the output'''
    name = 'VEPAnnotations'
    features = ['GeneID', 'FeatureID', 'GeneName', 'CCDS', 'Intron', 'Exon']
    VEP_feature_names = ['Gene', 'Feature', 'SYMBOL', 'CCDS', 'INTRON', 'EXON']
    consequence = True

    def process(self, res):
        for name, feature in zip(self.VEP_feature_names, self.features):
            if res[name] == '':
                del res[name]
                res[feature] = 'NA'
            else:
                res[feature] = res[name]
        return res

class VariantPosition(Annotation):
    name = 'VariantPosition'
    features = ['cDNApos','relcDNApos','CDSpos','relCDSpos','protPos','relProtPos']
    consequence = True

    def process(self, res):
        for num, raw_pos in enumerate(['cDNA_position', 'CDS_position', 'Protein_position']):
            if res[raw_pos] != '':
                pos, length = res[raw_pos].split('/')
                length = float(length)
                pos = pos.replace('?-','').replace('-?','').split('-')[0]
                if pos != '':
                    res[self.features[num*2]] = pos
                    res[self.features[num*2+1]] = ('%.3f' % (float(pos) / length)).rstrip('0').rstrip('.')
        return res


class oldVariantPosition(VariantPosition):
    name = 'oldVariantPosition'

    def process(self, res):
        for num, raw_pos in enumerate(['cDNA_position', 'CDS_position', 'Protein_position']):
            if res[raw_pos] != '':
                pos, length = res[raw_pos].split('/')
                length = float(length)
                pos = pos.replace('?-','').replace('-?','').split('-')[0]
                if pos != '':
                    res[self.features[num*2]] = pos
                    res[self.features[num*2+1]] = float(pos) / length
        return res

class Domain(Annotation):
    name = 'Domain'
    consequence = True
    domain_levels = ['ncoils', 'sigp', 'lcompl', 'ndomain', 'hmmpanther', 'other', 'UD']

    def process(self, res):
        if res['DOMAINS'] != '':
            domain_ind = 5
            for dfields in res['DOMAINS'].split('&'):
                try:
                    category,name = dfields.split(":", 1)
                except:
                    continue
                if category == "hmmpanther":
                    domain_ind = min(domain_ind, 4)
                elif ("_domain" in category) or ("_profile" in category):
                    domain_ind = min(domain_ind, 3)
                elif category == "Coiled-coils_(Ncoils)":
                    domain_ind = 0
                elif category == "Cleavage_site_(Signalp)":
                    domain_ind = min(domain_ind, 1)
                elif category == "Low_complexity_(Seg)":
                    domain_ind = min(domain_ind, 2)
            res['Domain'] = self.domain_levels[domain_ind]
        return res

class Dst2Splice(TabixAnnotation):
    name = 'Exon/Dst2Splice' # the loaded file is exon tabix
    features = ['Dst2Splice','Dst2SplType']
    path = '/reference/annotation_exons_v75.gtf.gz'
    consequence = True

    def process(self, res):
        '''This might be covered like in the normal tabix Annotations,
        unfortunately I currently cannot intergrate this, would
        need refactorization. Old code works though'''

        if res['FeatureID'].startswith('ENST'):
            minscore = None
            maxvalSPLICE = 20
            for iPos in range(res['Start'],res['End']+1):
                try:
                    for elem in self.tabix[0].fetch(res['Chrom'],
                                                 iPos-maxvalSPLICE-1,
                                                 iPos+maxvalSPLICE+1):
                        elem = elem.rstrip("\n").split('\t')
                        if res['FeatureID'] in elem[8]:
                            start,end = int(elem[3]),int(elem[4])
                            strand = elem[6] == "+"
                            if iPos < start:
                                val = (start - iPos,-1,True if strand else False)
                            elif iPos >= start and iPos <= end:
                                val =  min((iPos - start + 1,1,True if strand else False),(end - iPos + 1,1,False if strand else True))
                            else:
                                val = (iPos - end,-1,False if strand else True)
                            if (val[0] <= maxvalSPLICE) and ((minscore == None) or (val[0] < minscore[0])):
                                minscore = val
                                if minscore[0] == 1: break
                except: pass
            if not minscore is None:
                res["Dst2Splice"]=minscore[0]*minscore[1]
                res["Dst2SplType"]= "ACCEPTOR" if minscore[2] else "DONOR"
        return res

class MinDistTSX(TabixAnnotation):
    name = 'Transcript/MinDistTSX'
    features = ['minDistTSS', 'minDistTSE']
    path = '/reference/annotation_transcripts_v75.tsv.gz'
    consequence = True

    def process(self, res):

        minDistTSS,minDistTSE = None,None
        maxvalDIST = 10000000
        for iPos in set([res['Start'],res['End']]):
            for cmaxval in (maxvalDIST/10000,maxvalDIST/1000,maxvalDIST/100,maxvalDIST/10,maxvalDIST):
                try:
                    for elem in self.tabix[0].fetch(reference=res['Chrom'], start=max(0,iPos-cmaxval-1), end=iPos+cmaxval+1):
                        #sys.stderr.write("%s\n"%(elem))
                        elem = elem.rstrip("\n").split('\t')
                        start,end = int(elem[1]),int(elem[2])
                        if elem[3] == "-": start,end = end,start
                        val = min(max(start-iPos+1,iPos-start),maxvalDIST)
                        if minDistTSS == None or val < minDistTSS:
                            minDistTSS = val
                        val = min(max(end-iPos+1,iPos-end),maxvalDIST)
                        if minDistTSE == None or val < minDistTSE:
                            minDistTSE = val
                except: pass
                if minDistTSS != None and minDistTSE != None: break

        if minDistTSS != None:
            res["minDistTSS"] = str(minDistTSS)
        else:
            res["minDistTSS"] = str(maxvalDIST)
        if minDistTSE != None:
            res["minDistTSE"] = str(minDistTSE)
        else:
            res["minDistTSE"] = str(maxvalDIST)
        return res

class SIFT(Annotation):
    name = 'SIFT'
    features = ['SIFTcat', 'SIFTval']
    consequence = True

    def process(self, res):
        if res['SIFT'] != '':
            cat, val = res['SIFT'].split('(')
            res['SIFTval'] = val[:-1] # remove closing ')'
            res['SIFTcat'] = cat.split('_low_confidence')[0] # ignoring confidence for now
        return res

class PolyPhen(Annotation):
    name = 'PolyPhen'
    features = ['PolyPhenCat', 'PolyPhenVal']
    consequence = True

    def process(self, res):
        if res['PolyPhen'] != '':
            cat, val = res['PolyPhen'].split('(')
            res['PolyPhenVal'] = val[:-1] # remove closing ')'
            res['PolyPhenCat'] = cat
        return res

class PriPhCons(TabixAnnotation, ScoreHighest):
    name = 'priPhCons'
    path = '/phastCons/primates_nohuman.tsv.gz'

class MamPhCons(TabixAnnotation, ScoreHighest):
    name = 'mamPhCons'
    path = '/phastCons/placentalMammals_nohuman.tsv.gz'

class VerPhCons(TabixAnnotation, ScoreHighest):
    name = 'verPhCons'
    path = '/phastCons/vertebrates_nohuman.tsv.gz'

class PhastCons(TabixAnnotation, ScoreHighestPerFeature): # summarises the previous three in one Annotation
    name = 'PhastCons'
    features = ['priPhCons', 'mamPhCons', 'verPhCons']
    path = ''
    rangescore = True
    zerobased = True

class PriPhyloP(TabixAnnotation, ScoreHighest):
    name = 'priPhyloP'
    path = '/phyloP/primates_nohuman.tsv.gz'

class MamPhyloP(TabixAnnotation, ScoreHighest):
    name = 'mamPhyloP'
    path = '/phyloP/placentalMammals_nohuman.tsv.gz'

class VerPhyloP(TabixAnnotation, ScoreHighest):
    name = 'verPhyloP'
    path = '/phyloP/vertebrates_nohuman.tsv.gz'

class PhyloP(TabixAnnotation, ScoreHighestPerFeature): # summarises the previous three in one Annotation
    name = 'PhyloP'
    features = ['priPhyloP', 'mamPhyloP', 'verPhyloP']
    path = ''
    rangescore = True
    zerobased = True

class Bscores(TabixAnnotation, ScoreLowest):
    name = 'bStatistic'
    path = '/bscores/bscores.tsv.gz'
    multirange = True
    rangescore = True

class MirTargetScan(TabixAnnotation, ScoreHighest):
    name = 'targetScan'
    path = '/mirTargetScan/targetScanS.tsv.gz'
    multirange = True
    rangescore = True
    zerobased = True

class MirSVR(TabixAnnotation):
    name = 'mirSVR'
    path = '/mirSVR/mirSVR.tsv.gz'
    multirange = True
    rangescore = True
    features = ['mirSVR-Score','mirSVR-E','mirSVR-Aln']

    def _get_score(self, res):
        if len(self.score) == 1:
            res['mirSVR-Score'] = self.score[0][-1]
            res['mirSVR-E'] = self.score[0][-2]
            res['mirSVR-Aln'] = self.score[0][-3]
        elif len(self.score) > 1:
            helper = map(lambda x:(float(x[-1]),x[-3:]),self.score)
            helper.sort()
            res['mirSVR-Score'] = helper[0][-1][-1]
            res['mirSVR-E'] = helper[0][-1][-2]
            res['mirSVR-Aln'] = helper[0][-1][-3]
        return res

class ChromHMM(TabixAnnotation, ScoreAveragePerFeature):
    name = 'chromHMM_Epigenomes'
    path = '/chromHMM/roadMap_chromHMM.tsv.gz'
    features = ['cHmmTssA','cHmmTssAFlnk','cHmmTxFlnk','cHmmTx','cHmmTxWk',
                'cHmmEnhG','cHmmEnh','cHmmZnfRpts','cHmmHet','cHmmTssBiv',
                'cHmmBivFlnk','cHmmEnhBiv','cHmmReprPC','cHmmReprPCWk','cHmmQuies']
    rangescore = True

class ChromHMM_25state(TabixAnnotation, ScoreAveragePerFeature):
    name = 'chromHMM_25state'
    features = ['cHmm_E' + str(i) for i in range(1,26)]
    rangescore = True
    zerobased = True

class GERPelements(TabixAnnotation):
    name = 'GERPelements'
    path = '/gerp/gerp_elements.tsv.gz'
    rangescore = True
    features = ['GerpRS', 'GerpRSpval']

    def _get_score(self, res):
        if len(self.score) == 1:
            res['GerpRS'] = self.score[0][-2]
            res['GerpRSpval'] = self.score[0][-1]
        elif len(self.score) > 1:
            helper = map(lambda x:(float(x[-2]),x[-1],x[-2]),self.score)
            helper.sort()
            res['GerpRS'] = helper[0][-2]
            res['GerpRSpval'] = helper[-1][-1]
        return res

class GERPelementsBed(GERPelements):
    name = 'GERPelementsBed'
    path = ''
    zerobased = True

class GERPscores(TabixAnnotation):
    name = 'GERPscores'
    path = '/gerp/gerp_scores.tsv.gz'
    features = ['GerpN', 'GerpS']

    def _get_score(self, res):
        if len(self.score) == 1:
            res['GerpN'] = self.score[0][-2]
            res['GerpS'] = self.score[0][-1]
        elif len(self.score) > 1:
            helper = map(lambda x:(float(x[-2]),x[-2]),self.score)
            helper.sort()
            res['GerpN'] = helper[-1][-1]
            helper = map(lambda x:(float(x[-1]),x[-1]),self.score)
            helper.sort()
            res['GerpS'] = helper[-1][-1]
        return res

class GERPscoresBed(GERPscores):
    name = 'GERPscoresBed'
    path = ''
    rangescore = True
    zerobased = True

class Tfbs(TabixAnnotation):
    name = 'tfbs'
    path = '/tfbs/encodeAwgTfbsUniform.tsv.gz'
    features = ['TFBS', 'TFBSPeaks', 'TFBSPeaksMax']

    def _get_score(self, res):
        if len(self.score) > 0:
            for tfbsline in self.score:
                if ('TFBS' not in res) or (int(res['TFBS']) < int(tfbsline[-6])):
                    res['TFBS'] = tfbsline[-6]
                    res['TFBSPeaks'] = tfbsline[-4]
                    res['TFBSPeaksMax'] = tfbsline[-2]
        return res

class Motifs(TabixAnnotation):
    name = 'motifs'
    path = '/encodeMotifDB/motifs.tsv.gz'
    features = ['tOverlapMotifs','motifDist']

    def _get_score(self, res):
        if len(self.score) > 0:
            val1,val2 = 0,0
            failed = False
            for motif in self.score:
                sites = motif[-1].split(';')
                res['tOverlapMotifs'] = len(sites)
                ACGT = [0.0,0.0,0.0,0.0]
                for site in sites:
                    for ind,value in enumerate(site.split('_')[-4:]):
                        ACGT[ind] += float(value[1:])

                if res['Type'] in ['INS','DEL']:
                    val2 += max(ACGT)/len(sites)
                else:
                    if res['Ref'] == 'A': val1+=ACGT[0]/len(sites)
                    elif res['Ref'] == 'C': val1+=ACGT[1]/len(sites)
                    elif res['Ref'] == 'G': val1+=ACGT[2]/len(sites)
                    elif res['Ref'] == 'T': val1+=ACGT[3]/len(sites)
                    else: failed = True

                    if res['Alt'] == 'A': val2+=ACGT[0]/len(sites)
                    elif res['Alt'] == 'C': val2+=ACGT[1]/len(sites)
                    elif res['Alt'] == 'G': val2+=ACGT[2]/len(sites)
                    elif res['Alt'] == 'T': val2+=ACGT[3]/len(sites)
                    else: failed = True

            if not failed:
                res['motifDist']='%.2f'%(val1-val2)
        return res

# SEGWAY SCORING
segway_order = ['D','L0','L1','F0','F1','R0','R1','R2','R3','R4','R5','C0','C1','H3K9me1','TF0','TF1','TF2','TSS','GS','GE0','GE1','GE2','GM0','GM1','E/GM'][::-1]
segway_ranking = dict(zip(segway_order,range(len(segway_order))))
class Segway(TabixAnnotation):
    name = 'Segway'
    path = '/segway/segway.tsv.gz'
    rangescore = True

    def _get_score(self, res):
        if len(self.score) == 1:
            res[self.name] = self.score[0][-1]
        elif len(self.score) > 1:
            helper = [(segway_ranking[elem[-1]],elem[-1]) for elem in self.score]
            helper.sort()
            res[self.name] = helper[0][-1]
        return res

class EncodeH3K27Ac(EncodeAnnotation):
    name = 'EncH3K27Ac'
    path = '/H3K27Ac/tabix_files.lst'

class EncodeH3K4Me1(EncodeAnnotation):
    name = 'EncH3K4Me1'
    path = '/H3K4Me1/tabix_files.lst'

class EncodeH3K4Me3(EncodeAnnotation):
    name = 'EncH3K4Me3'
    path = '/H3K4Me3/tabix_files.lst'

class EncodeExp(EncodeAnnotation):
    name = 'EncExp'
    path = '/expression/tabix_files.lst'

class EncodeNucleosome(EncodeAnnotation):
    name = 'EncNucleo'
    path = '/nucleosomes/tabix_files.lst'

class EncodeOpenChrom(EncodeAnnotation):
    name = 'OpenChrom'
    path = '/openChromatin/tabix_files.lst'
    features = ['EncOCC','EncOCCombPVal','EncOCDNasePVal','EncOCFairePVal',
                'EncOCpolIIPVal','EncOCctcfPVal','EncOCmycPVal','EncOCDNaseSig',
                'EncOCFaireSig','EncOCpolIISig','EncOCctcfSig','EncOCmycSig']

    def _retrieve(self, res):
        #0: chrom   chr1    varchar(255)    values  Name of the chromosome
        #1: chromStart      10166   int(10) unsigned        range   Start position in chromosome
        #2: chromEnd        10376   int(10) unsigned        range   End position in chromosome
        #3: name    FAIREOnly_1     varchar(255)    values  Optional. Name given to a region (preferably unique). Use . if no name is assigned.
        #4: score   1000    int(10) unsigned        range   Optional. Indicates how dark the peak will be displayed in the browser (1-1000). If '0', the DCC will assign this based on signal value. Ideally average signalValue per base spread between 100-1000.
        #5: strand  .       char(1)         values  Optional. +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
        #6: thickStart      10166   int(10) unsigned        range   Start of where display should be thick (start codon)
        #7: thickEnd        10376   int(10) unsigned        range   End of where display should be thick (stop codon)
        #8: color   10027008        int(10) unsigned        range   RGB color value for peak
        #9: pValue  0.99    float unsigned  range   Fisher's Combined P-Value (-log 10)
        #10: dnaseSignal     0.0016  float unsigned  range   DNaseI HS Signal
        #11: dnasePvalue     0       float unsigned  range   DNaseI HS P-Value (-log 10)
        #12: faireSignal     0.0158  float unsigned  range   FAIRE Signal
        #13: fairePvalue     1.68    float unsigned  range   FAIRE P-Value (-log 10)
        #14: polIISignal     0.0265  float   range   Pol-II Signal
        #15: polIIPvalue     1.51    float   range   Pol-II P-Value (-log 10)
        #16: ctcfSignal      0.0084  float   range   CTCF Signal
        #17: ctcfPvalue      0       float   range   CTCF P-Value (-log 10)
        #18: cmycSignal      0.0336  float   range   c-Myc Signal
        #19: cmycPvalue      0       float   range   c-Myc P-Value (-log 10)
        #20: ocCode  3       int(10) unsigned        range   Open Chromatin (OC) Code

        maxscores = [[],[],[],[],[],[],[],[],[],[],[],[]] # [(20,locCode),(9,lcomb_pvalue),(11,ldnase_pvalue),(13,lfaire_pvalue),(15,lpolII_pvalue),(17,lctcf_pvalue),(19,lcmyc_pvalue),(10,ldnase_signal),(12,lfaire_signal),(14,lpolII_signal),(16,lctcf_signal),(18,lcmyc_signal)]
        newlist = []
        for encodeTabix in self.tabixlist:
            encodeTabix = get_range_from_tabix(encodeTabix,
                                                          res['Chrom'],
                                                          res['Start'],
                                                          res['End'],
                                                          rangescore=True,
                                                          continuous=self.continuous
                                                          )
            newlist.append(encodeTabix)
            encode_out = encodeTabix[6]
            if len(encode_out) > 0:
                #print encode_out
                for pos,ind in enumerate([20,9,11,13,15,17,19,10,12,14,16,18]):
                    if ind != 20:
                        value = maxL(map(lambda x:float(x[ind]),encode_out))
                    else:
                        value = minL(map(lambda x:float(x[ind]),encode_out))
                    if value != None: maxscores[pos].append(value)
        self.tabix = newlist
        self.score = map(maxL,maxscores)

    def _get_score(self, res):
        #ocCode,comb_pvalue,dnase_pvalue,faire_pvalue,polII_pvalue,ctcf_pvalue,cmyc_pvalue,dnase_signal,faire_signal,polII_signal,ctcf_signal,cmyc_signal = self.score
        if self.score[0] != None: res['EncOCC'] = self.score[0] # '%d'%ocCode
        if self.score[1] != None: res['EncOCCombPVal'] = self.score[1] # '%.3f'%comb_pvalue
        if self.score[2] != None: res['EncOCDNasePVal'] = self.score[2] # '%.3f'%dnase_pvalue
        if self.score[3] != None: res['EncOCFairePVal'] = self.score[3] # '%.3f'%faire_pvalue
        if self.score[4] != None: res['EncOCpolIIPVal'] = self.score[4] # '%.3f'%polII_pvalue
        if self.score[5] != None: res['EncOCctcfPVal'] = self.score[5] # '%.3f'%ctcf_pvalue
        if self.score[6] != None: res['EncOCmycPVal'] = self.score[6] # '%.3f'%cmyc_pvalue
        if self.score[7] != None: res['EncOCDNaseSig'] = self.score[7] # '%.2f'%dnase_signal
        if self.score[8] != None: res['EncOCFaireSig'] = self.score[8] # '%.2f'%faire_signal
        if self.score[9] != None: res['EncOCpolIISig'] = self.score[9] # '%.2f'%polII_signal
        if self.score[10] != None: res['EncOCctcfSig'] = self.score[10] # '%.2f'%ctcf_signal
        if self.score[11] != None: res['EncOCmycSig'] = self.score[11] # '%.2f'%cmyc_signal
        return res

class EncodeH3K4me1FC(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeH3K4me1FC' # FC = Fold Change (over control)
    features = ['EncodeH3K4me1-sum', 'EncodeH3K4me1-max']
    rangescore = True
    zerobased = True

class EncodeH3K4me2FC(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeH3K4Me2FC' # FC = Fold Change (over control)
    features = ['EncodeH3K4me2-sum', 'EncodeH3K4me2-max']
    rangescore = True
    zerobased = True

class EncodeH3K4me3FC(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeH3K4me3FC' # FC = Fold Change (over control)
    features = ['EncodeH3K4me3-sum', 'EncodeH3K4me3-max']
    rangescore = True
    zerobased = True

class EncodeH3K9acFC(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeH3K9acFC' # FC = Fold Change (over control)
    features = ['EncodeH3K9ac-sum', 'EncodeH3K9ac-max']
    rangescore = True
    zerobased = True

class EncodeH3K9me3FC(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeH3K9me3FC' # FC = Fold Change (over control)
    features = ['EncodeH3K9me3-sum', 'EncodeH3K9me3-max']
    rangescore = True
    zerobased = True

class EncodeH3K27acFC(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeH3K27acFC' # FC = Fold Change (over control)
    features = ['EncodeH3K27ac-sum', 'EncodeH3K27ac-max']
    rangescore = True
    zerobased = True

class EncodeH3K27me3FC(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeH3K27me3FC' # FC = Fold Change (over control)
    features = ['EncodeH3K27me3-sum', 'EncodeH3K27me3-max']
    rangescore = True
    zerobased = True

class EncodeH3K36me3FC(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeH3K36me3FC' # FC = Fold Change (over control)
    features = ['EncodeH3K36me3-sum', 'EncodeH3K36me3-max']
    rangescore = True
    zerobased = True

class EncodeH3K79me2FC(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeH3K79me2FC' # FC = Fold Change (over control)
    features = ['EncodeH3K79me2-sum', 'EncodeH3K79me2-max']
    rangescore = True
    zerobased = True

class EncodeH4K20me1FC(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeH4K20me1FC' # FC = Fold Change (over control)
    features = ['EncodeH4K20me1-sum', 'EncodeH4K20me1-max']
    rangescore = True
    zerobased = True

class EncodeH2AFZFC(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeH2AFZFC' # FC = Fold Change (over control)
    features = ['EncodeH2AFZ-sum', 'EncodeH2AFZ-max']
    rangescore = True
    zerobased = True

class EncodeDNase(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodeDNase'
    features = ['EncodeDNase-sum', 'EncodeDNase-max']
    rangescore = True
    zerobased = True

class EncodetotalRNA(TabixAnnotation, ScoreHighestPerFeature):
    name = 'EncodetotalRNA'
    features = ['EncodetotalRNA-sum', 'EncodetotalRNA-max']
    rangescore = True
    zerobased = True

class Grantham(AACombination):
    name = 'Grantham'
    path = '/grantham/grantham_matrix.tsv'
    consequence = True

class NeighboringMutations(TabixAnnotation):
    name = 'NeighboringMutations' # this annotation is replacement of the previous one (Dist2Mutation)
    features = ['Dist2Mutation']
    rangescore = True
    zerobased = True
    multirange = True

    def _get_score(self, res):
        if self.score:
            # distance between next mutation downstream and next mutation upstream
            # ignore mutations in the same position (or bridged ones)
            res['Dist2Mutation'] = int(self.score[-1][2]) - int(self.score[0][1]) - 1
        return res

class MutationDensity100(TabixAnnotation, ScoreHighestPerFeature):
    name = 'MutationDensity100'
    features = ['Freq100bp', 'Rare100bp', 'Sngl100bp']
    rangescore = True
    zerobased = True
    path = ''
    datatype = int

class MutationDensity1000(TabixAnnotation, ScoreHighestPerFeature):
    name = 'MutationDensity1000'
    features = ['Freq1000bp', 'Rare1000bp', 'Sngl1000bp']
    rangescore = True
    zerobased = True
    path = ''
    datatype = int

class MutationDensity10000(TabixAnnotation, ScoreHighestPerFeature):
    name = 'MutationDensity10000'
    features = ['Freq10000bp', 'Rare10000bp', 'Sngl10000bp']
    rangescore = True
    zerobased = True
    path = ''
    datatype = int

regElement_order = ['Promoter', 'TF binding site', 'Enhancer', 'Promoter Flanking Region', 'CTCF Binding Site', 'Open chromatin']
regElement_ranking = dict(zip(regElement_order,range(len(regElement_order))))
class RegulatoryBuiltElements(TabixAnnotation):
    name = 'EnsembleRegulatoryFeature'
    path = '/regulatoryBuilt/regulatory_features.bed.bgz'
    rangescore = True
    zerobased = True

    def _get_score(self, res):
        if len(self.score) == 1:
            res[self.name] = self.score[0][-1]
        elif len(self.score) > 1:
            helper = [(regElement_ranking[elem[-1]],elem[-1]) for elem in self.score]
            helper.sort()
            res[self.name] = helper[0][-1]
        return res

class DbscSNV(TabixAnnotation):
    name = 'dbscSNV'
    features = ['dbscSNV-ada_score', 'dbscSNV-rf_score']
    path = '/dbscSNV/dbscSNV1.1_GRCh37.txt.gz'

    def _get_score(self, res):
        for hit in self.score:
            if (hit[2] == res['Ref'] and hit[3] == res['Alt']) or (hit[3] == res['Ref'] and hit[2] == res['Alt']):
                res['dbscSNV-ada_score'] = hit[-2]
                res['dbscSNV-rf_score'] = hit[-1]
                break
        return res

class RemapOverlap(TabixAnnotation, ScoreHighestPerFeature):
    name = 'RemapOverlap'
    features = ['RemapOverlapTF', 'RemapOverlapCL']
    path = '/Remap/ReMap2_overlap_hg19.bg.gz'
    rangescore = True
    zerobased = True
    datatype = int

annotations = [
    Transversion(),
    Length(),
    Consequence(),
    GC_CpG(),
    oldGC_CpG(),
    MotifScores(),
    AminoAcids(),
    VEPAnnotations(),
    VariantPosition(),
    oldVariantPosition(),
    Domain(),
    Dst2Splice(),
    MinDistTSX(),
    SIFT(),
    PolyPhen(),
    PriPhCons(),
    MamPhCons(),
    VerPhCons(),
    PhastCons(),
    PriPhyloP(),
    MamPhyloP(),
    VerPhyloP(),
    PhyloP(),
    Bscores(),
    MirTargetScan(),
    MirSVR(),
    ChromHMM(),
    ChromHMM_25state(),
    GERPelements(),
    GERPelementsBed(),
    GERPscores(),
    GERPscoresBed(),
    Tfbs(),
    Motifs(),
    Segway(),
    EncodeH3K27Ac(),
    EncodeH3K4Me1(),
    EncodeH3K4Me3(),
    EncodeExp(),
    EncodeNucleosome(),
    EncodeOpenChrom(),
    EncodeH3K4me1FC(),
    EncodeH3K4me2FC(),
    EncodeH3K4me3FC(),
    EncodeH3K9acFC(),
    EncodeH3K9me3FC(),
    EncodeH3K27acFC(),
    EncodeH3K27me3FC(),
    EncodeH3K36me3FC(),
    EncodeH3K79me2FC(),
    EncodeH4K20me1FC(),
    EncodeH2AFZFC(),
    EncodeDNase(),
    EncodetotalRNA(),
    Grantham(),
    NeighboringMutations(),
    MutationDensity100(),
    MutationDensity1000(),
    MutationDensity10000(),
    RegulatoryBuiltElements(),
    DbscSNV(),
    RemapOverlap(),
    ]
