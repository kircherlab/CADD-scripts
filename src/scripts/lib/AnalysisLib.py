#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *08.05.2012
"""

import sys
import os
from collections import defaultdict
import pysam
import string
import mmap
import gzip
import traceback

verbose = False

table = string.maketrans('TGCABDHKMRSVWYtgcabdhkmrsvwy', 'ACGTVHDMKYSBWRACGTVHDMKYSBWR')  # COMPLEMENT DNA IGNORING CASE

# TRANSLATION OF EPO 6 PRIMATE SPECIES IDENTIFIERS
translations = {'Homo_sapiens': 'Hsap',
                'Pan_troglodytes': 'Ptro',
                'Gorilla_gorilla': 'Ggor',
                'Pongo_pygmaeus': 'Ppyg',
                'Macaca_mulatta': 'Mmul',
                'Callithrix_jacchus': 'Cjac',
                'Ancestral_sequences': 'Aseq'}


# COMPLEMENT DNA SEQUENCE
def compl(seq):
    global table
    return seq.translate(table)

# REVERSE COMPLEMENT DNA SEQUENCE


def revcompl(seq):
    global table
    seq = seq.translate(table)
    return seq[::-1]

# REVERSE DNA SEQUENCE


def reverse(seq):
    seq = list(seq)
    seq.reverse()
    return "".join(seq)

# SHORTEN SEQUENCE (USEFUL FOR ABBREVIATING LONG INDELS)


def shorten_seq(seq):
    if len(seq) > 12:
        sub_seq = seq.split(",")
        if len(sub_seq) == 1:
            return "%s..%d..%s" % (seq[:4], len(seq)-8, seq[-4:])
        else:
            res = []
            for seq in sub_seq:
                if len(seq) > 12:
                    res.append("%s..%d..%s" % (seq[:4], len(seq)-8, seq[-4:]))
                else:
                    res.append(seq)
            return ",".join(res)
    else:
        return seq

# SIMPLE FASTA PARSER


def read_fasta(filename):
    if filename == None:
        infile = sys.stdin
    else:
        if filename.endswith(".gz"):
            infile = gzip.open(filename)
        else:
            infile = open(filename)
    seq = ""
    header = ""
    for line in infile:
        line = line.strip()
        if (line[0] == ">"):
            if (header != "") and (seq != ""):
                yield header, seq
            header = line[1:]
            seq = ""
        else:
            seq += line
    if (header != ""):
        yield header, seq
    if filename != None:
        infile.close()
    raise StopIteration

# SIMPLE FASTQ PARSER


def read_fastq_file(filename):
    if filename == None:
        infile = sys.stdin
    elif filename.endswith(".gz"):
        infile = gzip.open(filename)
    else:
        infile = open(filename)
    name, seq, qual = "", "", ""
    count = 0
    for line in infile:
        if line.startswith("@") and ((count == 0) or (count == 3)):
            count = 1
            if len(seq) > 0:
                yield name, seq, qual
            name = line[1:].strip()
            seq = ""
            qual = ""
        elif (not line.startswith("+")) and (count == 1):
            seq += line.strip()
        elif line.startswith("+") and (count == 1):
            count = 2
        elif (count >= 2):
            qual += line.strip()
            count = 3
        else:
            if line.startswith('>'):
                for header, seq in read_fasta(filename):
                    yield header, seq, None
                raise StopIteration
            else:
                print("Unexpected line:", line.strip())
                raise StopIteration
    if len(seq) > 0:
        yield name, seq, qual
    if filename != None:
        infile.close()
    raise StopIteration


#################################################
###
# WORKING WITH MMAPPED FASTA GENOMES
###
#################################################

def read_genome_fastaindex(faifile):
    # READ FASTA INDEX TO MEMORY
    fastaindex = {}
    if os.path.exists(faifile):
        infile = open(faifile)
        for line in infile:
            fields = line.split()
            if len(fields) == 5:
                cname, length, start, line, cline = fields
                fastaindex[cname] = int(length), int(start), int(line), int(cline)
            else:
                sys.stderr.write('Error: Unexpected line in fasta index file: %s\n' % (line.strip()))
                sys.exit()
        infile.close()
    else:
        sys.stderr.write('Error: Fasta index file not available: %s\n' % faifile)
        sys.exit()
    return fastaindex


def get_DNA_file(chrom, start, end, fastafile, fastaindex):
    try:
        length, bstart, bline, cline = fastaindex[chrom]
        bases = ""
        for pos in xrange(start, end+1):
            if pos <= length and pos >= 1:
                hpos = pos - 1
                npos = (hpos // bline)*cline+(hpos % bline)
                fastafile.seek(bstart+npos)
                bases += fastafile.read(1)
        return bases
    except:
        return None


def get_DNA_iter(chrom, start, end, fastafile, fastaindex):
    try:
        length, bstart, bline, cline = fastaindex[chrom]
        for pos in xrange(start, end+1):
            if pos <= length and pos >= 1:
                hpos = pos - 1
                npos = (hpos // bline)*cline+(hpos % bline)
                fastafile.seek(bstart+npos)
                yield fastafile.read(1)
    except:
        pass
    raise StopIteration


def get_fasta_mmap(fastafile, length=0):
    cmap = None
    if os.path.exists(fastafile):
        # OPEN FASTA WITH MMAP
        f = open(fastafile, "r")
        cmap = mmap.mmap(f.fileno(), length=length, access=mmap.ACCESS_READ)
    else:
        sys.stderr.write('Error: Fasta file not available: %s\n' % fastafile)
        sys.exit()
    return cmap


def DNAiterator(chrom, start, fastamap, fastaindex):
    try:
        length, bstart, bline, cline = fastaindex[chrom]
        pos = start
        while (pos <= length and pos >= 1):
            hpos = pos - 1
            npos = (hpos // bline)*cline+(hpos % bline)
            fastamap.seek(bstart+npos)
            yield fastamap.read(1)
            pos += 1
    except:
        pass
    raise StopIteration


def get_DNA(chrom, start, end, fastamap, fastaindex):
    try:
        length, bstart, bline, cline = fastaindex[chrom]
        bases = ""
        for pos in xrange(start, end+1):
            if pos <= length and pos >= 1:
                hpos = pos - 1
                npos = (hpos // bline)*cline+(hpos % bline)
                fastamap.seek(bstart+npos)
                bases += fastamap.read(1)
        return bases
    except:
        return None


def count_GC_CpG(chrom, start, end, fastamap, fastaindex):
    try:
        length, bstart, bline, cline = fastaindex[chrom]
        CpG, GC = 0, 0
        count = 0
        lbase = ''
        for pos in xrange(start, end+1):
            if pos <= length and pos >= 1:
                hpos = pos - 1
                npos = (hpos // bline)*cline+(hpos % bline)
                fastamap.seek(bstart+npos)
                base = fastamap.read(1)
                count += 1
                if base in 'GC':
                    GC += 1
                if lbase == 'C' and base == 'G':
                    CpG += 1
                lbase = base
        if count > 0:
            return GC/float(count), CpG/(count*0.5)
        else:
            return None, None
    except:
        return None, None


#################################################
###
# WORKING WITH TABIX FILES
###
#################################################

#GERPelementsTabix = None
# if options.GERPelements != "" and (not os.path.exists(options.GERPelements) or not os.path.exists(options.GERPelements+".tbi")):
    #sys.stderr.write("GERPelements annotation: Require valid path to compressed tabix file and tabix index file.\n")
# elif os.path.exists(options.GERPelements) and os.path.exists(options.GERPelements+".tbi"):
    #sys.stderr.write('GERPelements annotation variants tabix file (%s)...\n'%(options.GERPelements))
    #GERPelementsTabix = pysam.Tabixfile(options.GERPelements,'r'),None,None,None,"GERPelements annotation"

#GERPscoresTabix = None
# if options.GERPscores != "" and (not os.path.exists(options.GERPscores) or not os.path.exists(options.GERPscores+".tbi")):
    #sys.stderr.write("GERPscores annotation: Require valid path to compressed tabix file and tabix index file.\n")
# elif os.path.exists(options.GERPscores) and os.path.exists(options.GERPscores+".tbi"):
    #sys.stderr.write('GERPscores annotation variants tabix file (%s)...\n'%(options.GERPscores))
    #GERPscoresTabix = pysam.Tabixfile(options.GERPscores,'r'),None,None,None,"GERPscores annotation"

#GERPscores,GERPelements = None,None
#GERPscoresTabix,GERPscores = get_from_tabix(GERPscoresTabix,chrom,pos)
#GERPelementsTabix,GERPelements = get_from_tabix(GERPelementsTabix,chrom,pos,rangescore=True)

# VCF FIELDS
fchr = 0
fpos = 1
fdbsnp = 2
fref_allele = 3
falt_allele = 4
fgeno_qual = 5
fflag = 6
finfo = 7
fformat = 8
fvalues = 9


def eval_syserror_frequencies(variation_data_entries, cref='.', alt='.'):
    AF = None
    for variation_data in variation_data_entries:
        alts = variation_data[falt_allele].split(',')
        for calt in alts:
            if (((alt == ".") or ((len(alt) == 1) and (calt == alt)) or ((len(alt) != 1) and (len(calt) == len(alt)))) and
                    ((cref == ".") or (variation_data[fref_allele] == cref))):  # IDENTIFY ENTRY THAT MATCHES THE DESCRIBED EXAMPLE (SNVs)/Length Indels
                for elem in variation_data[finfo].split(';'):
                    if elem.startswith('OBSERVED='):
                        counts = map(float, elem.split("=")[1].split('/'))
                        if counts[0] > counts[1]:
                            counts[1] = counts[0]
                        AF = counts[0]/counts[1]
                        break
    return AF


def eval_ESP_frequencies(variation_data_entries, cref='.', alt='.'):
    AF, AFR_AF, EUR_AF = None, None, None
    for variation_data in variation_data_entries:
        alts = variation_data[falt_allele].split(',')
        for calt in alts:
            if (((alt == ".") or ((len(alt) == 1) and (calt == alt)) or ((len(alt) != 1) and (len(calt) == len(alt)))) and
                    ((cref == ".") or (variation_data[fref_allele] == cref))):  # IDENTIFY ENTRY THAT MATCHES THE DESCRIBED EXAMPLE (SNVs)/Length Indels
                for elem in variation_data[finfo].split(';'):
                    if elem.startswith('TAC='):
                        counts = map(int, elem.split("=")[1].split(','))
                        for ind, calt in enumerate(alts):
                            if ((alt == ".") or ((len(alt) == 1) and (calt == alt)) or ((len(alt) != 1) and (len(calt) == len(alt)))):
                                if AF == None:
                                    AF = counts[ind]/float(sum(counts))
                                else:
                                    AF += counts[ind]/float(sum(counts))
                    elif elem.startswith('AA_AC='):
                        counts = map(int, elem.split("=")[1].split(','))
                        for ind, calt in enumerate(alts):
                            if ((alt == ".") or ((len(alt) == 1) and (calt == alt)) or ((len(alt) != 1) and (len(calt) == len(alt)))):
                                if AFR_AF == None:
                                    AFR_AF = counts[ind]/float(sum(counts))
                                else:
                                    AFR_AF += counts[ind]/float(sum(counts))
                    elif elem.startswith('EA_AC='):
                        counts = map(int, elem.split("=")[1].split(','))
                        for ind, calt in enumerate(alts):
                            if ((alt == ".") or ((len(alt) == 1) and (calt == alt)) or ((len(alt) != 1) and (len(calt) == len(alt)))):
                                if EUR_AF == None:
                                    EUR_AF = counts[ind]/float(sum(counts))
                                else:
                                    EUR_AF += counts[ind]/float(sum(counts))
    return AF, AFR_AF, EUR_AF


def eval_1000G_frequencies(variation_data_entries, cref='.', alt='.', simple=False):
    if not simple:
        return eval_dbSNP_frequencies(variation_data_entries, cref, alt)
    else:
        duplicates = set()
        AF1000g = None
        for variation_data in variation_data_entries:
            alts = variation_data[falt_allele].split(',')
            for cind, calt in enumerate(alts):
                if (((alt == ".") or ((len(alt) == 1) and (calt == alt)) or ((len(alt) != 1) and (len(calt) == len(alt)))) and
                        ((cref == ".") or (variation_data[fref_allele] == cref))):  # IDENTIFY ENTRY THAT MATCHES THE DESCRIBED EXAMPLE (SNVs)/Length Indels
                    if (variation_data[falt_allele], variation_data[fref_allele]) in duplicates:
                        continue
                    duplicates.add((variation_data[falt_allele], variation_data[fref_allele]))
                    for elem in variation_data[finfo].split(';'):
                        if elem.startswith('AF='):
                            if AF1000g == None:
                                AF1000g = float(elem[3:].split(",")[cind])
                            else:
                                AF1000g += float(elem[3:].split(",")[cind])
        return AF1000g


def eval_dbSNP_frequencies(variation_data_entries, cref='.', alt='.'):
    duplicates = set()
    AF1000g, ASN_AF, AMR_AF, AFR_AF, EUR_AF, isLowCov = None, None, None, None, None, None
    for variation_data in variation_data_entries:
        alts = variation_data[falt_allele].split(',')
        for calt in alts:
            if (((alt == ".") or ((len(alt) == 1) and (calt == alt)) or ((len(alt) != 1) and (len(calt) == len(alt)))) and
                    ((cref == ".") or (variation_data[fref_allele] == cref))):  # IDENTIFY ENTRY THAT MATCHES THE DESCRIBED EXAMPLE (SNVs)/Length Indels
                if (variation_data[falt_allele], variation_data[fref_allele]) in duplicates:
                    continue
                duplicates.add((variation_data[falt_allele], variation_data[fref_allele]))
                isLowCov = False
                for elem in variation_data[finfo].split(';'):
                    if elem.startswith('AF='):
                        if AF1000g == None:
                            AF1000g = float(elem[3:])
                        else:
                            AF1000g += float(elem[3:])
                    elif elem.startswith('ASN_AF='):
                        if ASN_AF == None:
                            ASN_AF = float(elem[7:])
                        else:
                            ASN_AF += float(elem[7:])
                    elif elem.startswith('AMR_AF='):
                        if AMR_AF == None:
                            AMR_AF = float(elem[7:])
                        else:
                            AMR_AF += float(elem[7:])
                    elif elem.startswith('AFR_AF='):
                        if AFR_AF == None:
                            AFR_AF = float(elem[7:])
                        else:
                            AFR_AF += float(elem[7:])
                    elif elem.startswith('EUR_AF='):
                        if EUR_AF == None:
                            EUR_AF = float(elem[7:])
                        else:
                            EUR_AF += float(elem[7:])
                    elif elem.startswith("SNPSOURCE=") and "LOWCOV" in elem:
                        isLowCov = True
    return AF1000g, ASN_AF, AMR_AF, AFR_AF, EUR_AF, isLowCov


def get_from_tabix(tabixobj, chrom, pos, rangescore=False, chromcol=0, startcol=1, endcol=2, continuous=False, multirange=False, zerobased=False):
    if tabixobj == None:
        return tabixobj, []
    filestream, tabix_iter, tabix_next, tabix_chrom, old_start, old_end, old_res, name = tabixobj
    if (old_start == pos) and (tabix_chrom == chrom):
        # Return the same results again
        return tabixobj
    res = []
    if not rangescore:
        if (chrom != tabix_chrom) or ((chrom == tabix_chrom) and ((tabix_next[0] == None) or (pos < tabix_next[0] and not continuous) or (pos > tabix_next[0]+1000))):
            if verbose:
                sys.stderr.write('(Re-)Initiating iterator for %s tabix file...\n' % (name))
            tabix_chrom = chrom
            try:
                tabix_iter = filestream.fetch(reference=chrom, start=pos-1)
                fields = tabix_iter.next().rstrip('\n').split('\t')
                tabix_next = int(fields[startcol])+zerobased, fields
            except:
                #sys.stderr.write("Tabix exception for %s of %s %d\n"%(name,chrom,pos))
                tabix_next = None, None
        while (tabix_next[0] != None) and (pos >= tabix_next[0]):
            current = tabix_next
            try:
                fields = tabix_iter.next().rstrip('\n').split('\t')
                tabix_next = int(fields[startcol])+zerobased, fields
            except:
                if verbose:
                    sys.stderr.write("Tabix exception for %s of %s %d\n" % (name, chrom, pos))
                tabix_next = None, None
            if pos == current[0]:
                res.append(current[1])
    else:
        # Essentially impossible to save seeks operations if multiple intervals can overlap a position
        if (chrom != tabix_chrom) or multirange or ((chrom == tabix_chrom) and ((tabix_next[0] == None) or (pos < tabix_next[0] and not continuous) or (pos > int(tabix_next[1][endcol])+10000+zerobased))):
            if verbose:
                sys.stderr.write('(Re-)Initiating iterator for %s tabix file...\n' % (name))
            tabix_chrom = chrom
            try:
                tabix_iter = filestream.fetch(reference=chrom, start=pos-1)
                fields = tabix_iter.next().rstrip('\n').split('\t')
                tabix_next = int(fields[startcol])+zerobased, fields
            except:
                if verbose:
                    sys.stderr.write("Tabix exception for %s of %s %d\n" % (name, chrom, pos))
                tabix_next = None, None
        if tabix_next[0] != None and (pos >= tabix_next[0]) and (pos <= int(tabix_next[1][endcol])+zerobased):
            res.append(tabix_next[1])
        else:
            while (tabix_next[0] != None) and (pos >= int(tabix_next[1][endcol])+zerobased):
                try:
                    fields = tabix_iter.next().rstrip('\n').split('\t')
                    tabix_next = int(fields[startcol])+zerobased, fields
                except:
                    if verbose:
                        sys.stderr.write("Tabix exception for %s of %s %d\n" % (name, chrom, pos))
                    tabix_next = None, None
            if (tabix_next[0] != None) and (pos >= tabix_next[0]) and (pos < int(tabix_next[1][endcol])+zerobased):
                res.append(tabix_next[1])
                tabixobj = filestream, tabix_iter, tabix_next, tabix_chrom, pos, None, res, name
    return tabixobj


def get_range_from_tabix(tabixobj, chrom, start, end, rangescore=False, chromcol=0, startcol=1, endcol=2, continuous=False, multirange=False, zerobased=False):
    if tabixobj == None:
        return tabixobj, []
    filestream, tabix_iter, tabix_next, tabix_chrom, old_start, old_end, old_res, name = tabixobj
    if (old_start == start) and (old_end == end) and (tabix_chrom == chrom):
        # Return the same results again
        return tabixobj
    res = []
    pos = start
    if not rangescore:
        if (chrom != tabix_chrom) or ((chrom == tabix_chrom) and ((tabix_next[0] == None) or (pos < tabix_next[0] and not continuous) or (pos > tabix_next[0]+1000))):
            if verbose:
                sys.stderr.write('(Re-)Initiating iterator for %s tabix file...\n' % (name))
            tabix_chrom = chrom
            try:
                tabix_iter = filestream.fetch(reference=chrom, start=pos-1, end=end)
                for elem in tabix_iter:
                    fields = elem.rstrip('\n').split('\t')
                    res.append(fields)
                if continuous:
                    tabix_iter = filestream.fetch(reference=chrom, start=end)
                    fields = tabix_iter.next().rstrip('\n').split('\t')
                    tabix_next = int(fields[startcol])+zerobased, fields
                else:
                    tabix_next = None, None
            except:
                if verbose:
                    sys.stderr.write("Tabix exception for %s of %s %d\n" % (name, chrom, pos))
                tabix_next = None, None
        else:
            for pos in range(start, end+1):
                while (tabix_next[0] != None) and (pos >= tabix_next[0]):
                    current = tabix_next
                    try:
                        fields = tabix_iter.next().rstrip('\n').split('\t')
                        tabix_next = int(fields[startcol])+zerobased, fields
                    except:
                        if verbose:
                            sys.stderr.write("Tabix exception for %s of %s %d\n" % (name, chrom, pos))
                        tabix_next = None, None
                    if pos == current[0]:
                        res.append(current[1])
    else:
        # Essentially impossible to save seeks operations if multiple intervals can overlap a position
        if (chrom != tabix_chrom) or multirange or ((chrom == tabix_chrom) and ((tabix_next[0] == None) or (pos < tabix_next[0] and not continuous) or (pos > int(tabix_next[1][endcol])+10000+zerobased))):
            if verbose:
                sys.stderr.write('(Re-)Initiating iterator for %s tabix file...\n' % (name))
            tabix_chrom = chrom
            try:
                tabix_iter = filestream.fetch(reference=chrom, start=pos-1, end=end)
                for elem in tabix_iter:
                    fields = elem.rstrip('\n').split('\t')
                    res.append(fields)
                if continuous and not multirange:
                    tabix_iter = filestream.fetch(reference=chrom, start=end)
                    fields = tabix_iter.next().rstrip('\n').split('\t')
                    tabix_next = int(fields[startcol])+zerobased, fields
                else:
                    tabix_next = None, None
            except:
                if verbose:
                    sys.stderr.write("Tabix exception for %s of %s %d\n" % (name, chrom, pos))
                tabix_next = None, None
        else:
            for pos in range(start, end+1):
                if tabix_next[0] != None and (pos >= tabix_next[0]) and (pos < int(tabix_next[1][endcol])+zerobased):
                    res.append(tabix_next[1])
                else:
                    while (tabix_next[0] != None) and (pos >= int(tabix_next[1][endcol])+zerobased):
                        try:
                            fields = tabix_iter.next().rstrip('\n').split('\t')
                            tabix_next = int(fields[startcol])+zerobased, fields
                        except:
                            if verbose:
                                sys.stderr.write("Tabix exception for %s of %s %d\n" % (name, chrom, pos))
                            tabix_next = None, None
                    if (tabix_next[0] != None) and (pos >= tabix_next[0]) and (pos < int(tabix_next[1][endcol])+zerobased):
                        res.append(tabix_next[1])
    tabixobj = filestream, tabix_iter, tabix_next, tabix_chrom, start, end, res, name
    return tabixobj


#################################################
###
# WORKING WITH EPO ALIGNMENTS
###
#################################################

# EXAMPLE

#parser = OptionParser()
#parser.add_option("-a", "--ancestor_path", dest="ancestor_path", help="PATH of the EMF and EMF index files (def epo_6_primate_v66/split_mod)",default="epo_6_primate_v66/split_mod")
#parser.add_option("-r", "--reference", dest="reference", help="Reference species (def Hsap)",default="Hsap")
#(options, args) = parser.parse_args()

# OUTPUT ERROR MESSAGES AND EXIT IF THERE ARE MISSING FILES
# if not os.path.isdir(options.ancestor_path) or not os.path.exists(options.ancestor_path+"/hsa_emf.index") or not os.path.exists(options.ancestor_path+"/ptr_emf.index"):
    #sys.stderr.write("Require valid path to EMF and EMF index files.\n")
    # sys.exit()

# MAKE AN INDEX OF HUMAN EMF FILES
#hsa_index = emf_init_index(options.ancestor_path+"/hsa_emf.index")
# MAKE AN INDEX OF CHIMP EMF FILES
#ptr_index = emf_init_index(options.ancestor_path+"/ptr_emf.index")
# CREATE DATASTRUCTURE DESCRIBING CURRENT DATA BLOCK
# current = { "data":None, "filename":None, "species":None,
    # "chrom":None, "pos":None,
    # "block": None, "bline":None, "bstart":None, "bend":None,
    # "bseqs": None, "bpos": None, "bstrands": None, "bgaps": None,
    # "bancestors": None, "bref":None, "btype":None }
#emf_obj = (current,hsa_index,ptr_index,options.ancestor_path)

#if options.reference == "Hsap": CAnc = "SELF-Ptro"
# else: CAnc = "SELF-Hsap"

#hp = format_ancestor(emf_get_base_at_position(emf_obj,options.reference,chrom,pos,CAnc,True))
#hg = format_ancestor(emf_get_base_at_position(emf_obj,options.reference,chrom,pos,"SELF-Ggor",True))
#ho = format_ancestor(emf_get_base_at_position(emf_obj,options.reference,chrom,pos,"SELF-Pabe",True))
#nref = format_ancestor(emf_get_base_at_position(emf_obj,options.reference,chrom,pos+1,"SELF",True))
#nhp = format_ancestor(emf_get_base_at_position(emf_obj,options.reference,chrom,pos+1,CAnc,True))

# RETURN A DICTIONARY (CHROM) OF DICTIONARIES (POS) WITH INFORMATION ABOUT FILE, BYTES AND LINES
def emf_init_index(filename):
    from bx.intervals.intersection import Intersecter, Interval
    if os.path.isfile(filename):
        infile = open(filename)
        res = {}
        for line in infile:
            if line.startswith("#"):
                continue
            fields = line.split()  # Chr\tStart\tEnd\tFile\tByte\tLine
            if len(fields) == 6:
                chrom, start, end, cfile, bytes, lines = fields
                start, end = int(start)+1, int(end)+1  # MAKE COORDINATES ONE BASED HALF-OPEN
                if chrom not in res:
                    res[chrom] = {}
                    res[chrom]["coords"] = Intersecter()
                res[chrom]["coords"].add_interval(Interval(start, end))
                res[chrom][(start, end)] = cfile, int(bytes), int(lines)
            else:
                sys.stderr.write("Unexpected line [%s]: %s" % (filename, line))
        infile.close()
        return res
    return None

   # RETURN EMF BLOCK TYPE


def emf_get_type(emf_obj):
    current = emf_obj[0]
    if current["btype"] == None:
        return "NA"
    else:
        return current["btype"]


def emf_get_base_at_position(emf_obj, species, chrom, pos, outspecies="Hsap-Ptro", silent=False):
    current = emf_obj[0]
    hsa_index = emf_obj[1]
    ptr_index = emf_obj[2]
    ancestor_path = emf_obj[3]

    global translations
    if species == "Hsap":
        cindex = hsa_index
    elif species == "Ptro":
        cindex = ptr_index
    else:
        return None

    if chrom.startswith('chr'):
        chrom = chrom[3:]
    if chrom.startswith('gl'):
        chrom = chrom.upper()
    if chrom == "M":
        chrom = "MT"
    if chrom == "23":
        chrom = "X"
    result = []

    if current['species'] != species or current['chrom'] != chrom or (pos < current['bstart']) or (pos >= current['bend']):
        if chrom in cindex:
            res = cindex[chrom]["coords"].find(pos, pos+1)
            if verbose:
                print(chrom, pos, res)
            if len(res) == 1:
                rstart = res[0].start
                rend = res[0].end
                filename, bytes, loffset = cindex[chrom][(rstart, rend)]
                # bytes,loffset=int(bytes),int(loffset)
                if current['filename'] != filename:
                    if verbose:
                        sys.stderr.write("Reading file: %s\n" % (ancestor_path+"/"+filename))
                    if os.path.isfile(ancestor_path+"/"+filename) and filename.endswith('.gz'):
                        f = gzip.open(ancestor_path+"/"+filename, 'rb')
                        del current['data']
                        current['data'] = f.read().splitlines()
                        f.close()
                        current['filename'] = filename
                    else:
                        sys.stderr.write("Error: EMF input file (%s) is not available\n" % (filename))
                current['species'] = species
                current['chrom'] = chrom
                current['btype'] = None
                start_data = False
                is_compl = False
                found_pos = False
                current['block'] = []
                current['bline'] = None
                current['bstart'] = rstart
                current['bend'] = rend
                current['bseqs'] = []
                current['bpos'] = []
                current['bstrands'] = []
                current['bgaps'] = []
                current['bancestors'] = {}
                current['bref'] = None
                if verbose:
                    sys.stderr.write("Initiating block (%s %d-%d)...\n" % (chrom, rstart, rend))
                for line in current['data'][loffset:]:
                    if line.startswith('#'):
                        continue  # HEADER
                    elif line.startswith('SCORE'):
                        continue
                    elif line.startswith('//'):  # END OF ENTRY
                        break
                        #sys.stderr.write("Unexpected error: hit block end...\n")
                    elif line.startswith('SEQ'):
                        fields = line.split()
                        if fields[2] in translations:
                            fields[2] = translations[fields[2]]
                        # species, chromosome, start, end, strand, length
                        cid = (fields[1][:1].upper()+fields[1].split('_')[1][:3], fields[2], int(fields[3])-1,
                               int(fields[4]), ('-' if (fields[5][0] == '-') else '+'), int(fields[6].split('=')[1][:-1]))
                        if species == cid[0] and chrom == cid[1] and rstart == cid[2]+1 and rend == cid[3]+1:
                            current['bref'] = len(current['bseqs'])
                            current['bline'] = 0
                            is_compl = (fields[5][0] == '-')
                        current['bstrands'].append(fields[5][0] != '-')
                        current['bgaps'].append(0)
                        if current['bstrands'][-1]:
                            current['bpos'].append(int(fields[3])-1)
                        else:
                            current['bpos'].append(int(fields[4])+1)
                        current['bseqs'].append(cid)

                    elif line.startswith('TREE'):
                        refname = "%s_%s_%d_%d[%s]" % (current['bseqs'][current['bref']][0], current['bseqs'][current['bref']][1], current['bseqs']
                                                       [current['bref']][2]+1, current['bseqs'][current['bref']][3], current['bseqs'][current['bref']][4])
                        tree = line
                        keller = []
                        hname = ''
                        hpos = 0
                        while hpos < len(tree):
                            elem = tree[hpos]
                            if elem == '(':
                                if hname != '':
                                    keller.append(hname)
                                    hname = ''
                                keller.append(elem)
                            elif elem == ',':
                                if hname != '':
                                    keller.append(hname)
                                    hname = ''
                            elif elem == ')':
                                if hname != '':
                                    keller.append(hname)
                                    hname = ''
                                ancestor = ''
                                hpos += 1
                                while (hpos < len(tree)) and (tree[hpos] != ',') and (tree[hpos] != ')'):
                                    elem = tree[hpos]
                                    ancestor += elem
                                    hpos += 1
                                ancestor = ancestor.split(':')[0]
                                if (hpos < len(tree)):
                                    hpos -= 1
                                # print '--------------'
                                # print "K-I",keller
                                # print 'Ancestor',ancestor
                                value = keller.pop()
                                reorder = []
                                while value != '(':
                                    if (value.split(':')[0] == refname):
                                        value = '!'+value.replace(species, "SELF")
                                    elif value.startswith('SELF'):
                                        value = '#'+value
                                    elif value.startswith(species):
                                        value = '$'+value
                                    else:
                                        value = '-'+value
                                    reorder.append(value.split(':')[0])
                                    value = keller.pop()
                                reorder.sort()
                                val1 = reorder[0].split('_')[0][1:].split('-')[0]
                                val2 = reorder[-1].split('_')[0][1:].split('-')[0]
                                if val1 == val2:
                                    save = val1+'_'+'_'.join(ancestor.split('_')[1:][:-2])
                                else:
                                    save = val1+'-'+val2+'_'+'_'.join(ancestor.split('_')[1:][:-2])
                                keller.append(save)
                                current['bancestors']['_'.join(save.split('_')[1:])] = save.split('_')[0]
                                # print "K-O",keller
                            else:
                                hname += elem
                            hpos += 1
                    elif line.startswith('DATA'):
                        start_data = True
                        current['btype'] = ''.join(map(lambda x: x[0][0] if not x[0].startswith(
                            'Pabe') else 'O', filter(lambda x: not x[0].startswith('Aseq'), current['bseqs'])))
                        # print current['bancestors']
                        # print current['bseqs']
                        for hpos, elem in enumerate(current['bseqs']):
                            if elem[1] in current['bancestors']:
                                helper = list(elem)
                                helper[0] = current['bancestors'][elem[1]]
                                current['bseqs'][hpos] = tuple(helper)
                        # print current['bseqs']
                    elif start_data and len(line) >= len(current['bseqs']):
                        line = "".join(line.split(' '))  # FORMAT SPECIFICATION ALLOWS FOR SPACES
                        if is_compl:
                            line = compl(line)
                        current['block'].append(line)
                        # COUNT BASES:
                        if not found_pos:
                            for hpos, cid in enumerate(current['bseqs']):
                                if line[hpos] == "-":
                                    current['bgaps'][hpos] += 1
                                elif line[hpos] == "~":
                                    continue
                                else:
                                    current['bgaps'][hpos] = 0
                                    if current['bstrands'][hpos]:
                                        current['bpos'][hpos] += 1
                                    else:
                                        current['bpos'][hpos] -= 1
                            if current['bpos'][current['bref']] == pos:
                                found_pos = True
                                current['bline'] = len(current['block'])-1
                                #hhcount = 0
                        # if found_pos and hhcount < 20:
                            # print line.strip()
                            #hhcount += 1
            elif len(res) == 0:
                if not silent:
                    sys.stderr.write("Error: Coordinate (%s:%d) not available from EMF alignments\n" % (chrom, pos))
            else:
                sys.stderr.write("Error: Coordinate (%s:%d) seems present in several EMF alignments\n" % (chrom, pos))
        else:
            if not silent:
                sys.stderr.write("Error: Chromosome (%s) not available from EMF alignments\n" % (chrom))

    if (len(result) == 0) and (current['species'] == species) and (current['chrom'] == chrom) and (pos >= current['bstart']) and (pos < current['bend']):
        # GO TO POSITION
        # FORWARD
        if (current['bstrands'][current['bref']] and current['bpos'][current['bref']] <= pos) or (not current['bstrands'][current['bref']] and current['bpos'][current['bref']] >= pos):
            while (current['bline'] < (len(current['block'])-1)) and (current['bpos'][current['bref']] != pos):
                current['bline'] += 1
                line = current['block'][current['bline']]
                for hpos, cid in enumerate(current['bseqs']):
                    if line[hpos] == "-":
                        current['bgaps'][hpos] += 1
                    elif line[hpos] == "~":
                        continue
                    else:
                        current['bgaps'][hpos] = 0
                        if current['bstrands'][hpos]:
                            current['bpos'][hpos] += 1
                        else:
                            current['bpos'][hpos] -= 1
        # REVERSE
        else:
            while (current['bline'] > 0) and (current['bpos'][current['bref']] != pos):
                current['bline'] -= 1
                line = current['block'][current['bline']]
                for hpos, cid in enumerate(current['bseqs']):
                    if line[hpos] == "-":
                        current['bgaps'][hpos] -= 1
                    elif line[hpos] == "~":
                        continue
                    else:
                        current['bgaps'][hpos] = 0
                        if current['bstrands'][hpos]:
                            current['bpos'][hpos] -= 1
                        else:
                            current['bpos'][hpos] += 1
        # ARE AT THE POSITION
        if current['bpos'][current['bref']] == pos:
            # SAVE CURRENT POSITION
            helper = (list(current['bgaps']), list(current['bpos']), current['bline'])
            if current['bstrands'][current['bref']]:  # PLUS STRAND
                # WHILE WE ARE AT THE SAME REFERENCE POSITION WE COLLECT ALL BASES FROM OTHERS
                while (current['bline'] < len(current['block'])) and (current['bpos'][current['bref']] == pos):
                    line = current['block'][current['bline']]
                    if outspecies == "SELF":
                        if len(result) == 0:
                            result.append(line[current['bref']])
                        else:
                            result[0] += line[current['bref']]
                    elif outspecies == "ALL":
                        if len(result) == 0:
                            for hpos, cid in enumerate(current['bseqs']):
                                if not cid[0].startswith('Aseq') and not cid[0].startswith('SELF-') and not cid[1].startswith('Ancestor_'):
                                    result.append(line[hpos])
                        else:
                            ind = 0
                            for hpos, cid in enumerate(current['bseqs']):
                                if not cid[0].startswith('Aseq') and not cid[0].startswith('SELF-') and not cid[1].startswith('Ancestor_'):
                                    if line[hpos] != "-":
                                        result[ind] += line[hpos]
                                    ind += 1
                    else:
                        if len(result) == 0:
                            for hpos, cid in enumerate(current['bseqs']):
                                if cid[0] == outspecies:
                                    result.append(line[hpos])
                        else:
                            ind = 0
                            for hpos, cid in enumerate(current['bseqs']):
                                if cid[0] == outspecies:
                                    if line[hpos] != "-":
                                        result[ind] += line[hpos]
                                    ind += 1
                    current['bline'] += 1
                    if current['bline'] < len(current['block']):
                        line = current['block'][current['bline']]
                        for hpos, cid in enumerate(current['bseqs']):
                            if line[hpos] == "-":
                                current['bgaps'][hpos] += 1
                            elif line[hpos] == "~":
                                continue
                            else:
                                current['bgaps'][hpos] = 0
                                if current['bstrands'][hpos]:
                                    current['bpos'][hpos] += 1
                                else:
                                    current['bpos'][hpos] -= 1
            else:  # MINUS STRAND
                # WHILE WE ARE AT THE SAME REFERENCE POSITION WE COLLECT ALL BASES FROM OTHERS
                while (current['bline'] >= 0) and (current['bpos'][current['bref']] == pos):
                    line = current['block'][current['bline']]
                    if outspecies == "SELF":
                        if len(result) == 0:
                            result.append(line[current['bref']])
                        else:
                            result[0] += line[current['bref']]
                    elif outspecies == "ALL":
                        if len(result) == 0:
                            for hpos, cid in enumerate(current['bseqs']):
                                if not cid[0].startswith('Aseq') and not cid[0].startswith('SELF-') and not cid[1].startswith('Ancestor_'):
                                    result.append(line[hpos])
                        else:
                            ind = 0
                            for hpos, cid in enumerate(current['bseqs']):
                                if not cid[0].startswith('Aseq') and not cid[0].startswith('SELF-') and not cid[1].startswith('Ancestor_'):
                                    if line[hpos] != "-":
                                        result[ind] += line[hpos]
                                    ind += 1
                    else:
                        if len(result) == 0:
                            for hpos, cid in enumerate(current['bseqs']):
                                if cid[0] == outspecies:
                                    result.append(line[hpos])
                        else:
                            ind = 0
                            for hpos, cid in enumerate(current['bseqs']):
                                if cid[0] == outspecies:
                                    if line[hpos] != "-":
                                        result[ind] += line[hpos]
                                    ind += 1
                    current['bline'] -= 1
                    if (current['bline'] >= 0):
                        line = current['block'][current['bline']]
                        for hpos, cid in enumerate(current['bseqs']):
                            if line[hpos] == "-":
                                current['bgaps'][hpos] -= 1
                            elif line[hpos] == "~":
                                continue
                            else:
                                current['bgaps'][hpos] = 0
                                if current['bstrands'][hpos]:
                                    current['bpos'][hpos] -= 1
                                else:
                                    current['bpos'][hpos] += 1
            # SET BACK POSITION TO FIRST HIT
            current['bgaps'], current['bpos'], current['bline'] = helper
    return result


def format_ancestor(hlist, last=None, unique=True):
    if last != None:
        if len(hlist) == 0 and last == "N/A":
            return "N/A"
        elif len(hlist) == 0:
            return last
        elif len(hlist) == 1 and last == "N/A":
            return hlist[0].upper()
        elif len(hlist) == 1:
            return last+hlist[0].upper()
        elif len(hlist) > 1 and last == "N/A":
            if unique:
                return ",".join(list(set(map(lambda x: x.upper(), hlist))))
            else:
                return ",".join(map(lambda x: x.upper(), hlist))
        else:
            seqs = last.split(',')
            if len(seqs) == hlist:
                hlist = map(lambda x: x.upper(), hlist)
            else:
                hlist = list(set(map(lambda x: x.upper(), hlist)))
                if len(seqs) != hlist:
                    return "N/A"
            for elem, ind in seqs:
                hlist[ind] = elem+hlist[ind]
            return ",".join(hlist)
    else:
        if len(hlist) == 0:
            return "N/A"
        elif len(hlist) == 1:
            return hlist[0].upper()
        else:
            if unique:
                return ",".join(list(set(map(lambda x: x.upper(), hlist))))
            else:
                return ",".join(map(lambda x: x.upper(), hlist))
