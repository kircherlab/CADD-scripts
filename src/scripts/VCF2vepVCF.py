#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *20.11.2013
"""

import sys, os
ROOTdir = os.path.dirname(os.path.realpath(__file__))+"/../"
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
ROOTanno = (ROOTdir+"/annotations").replace("//",'/')

import pysam
from optparse import OptionParser

parser = OptionParser("%prog [options]")
parser.add_option("-s","--SNVs", dest="SNVs", help="Path to SNV output file (def STDOUT)",default=None)
parser.add_option("-i","--InDels", dest="InDels", help="Path to InDel output file (def STDOUT)",default=None)
parser.add_option("-r","--reference", dest="reference", help="Reference genome to perform reference check (fasta index needs to be available) (def %s/reference/reference.fa)"%(ROOTanno),default='%s/reference/reference.fa'%(ROOTanno))
parser.add_option("--refCheck", dest="refCheck", help="Perform reference check (def off)",default=False,action="store_true")
(options, args) = parser.parse_args()

reference = None
if options.refCheck:
  if not os.path.exists(options.reference) or not os.path.exists(options.reference+".fai"): 
    options.refCheck = False
    sys.stderr.write("Warning: Reference not available, can not perform reference check.\n")
  reference = pysam.Fastafile(options.reference)

def is_nucleotide(seq):
  for base in seq.upper():
    if base not in 'ACGT': return False
  return True

def sharedPrefix(s1,s2):
  minLength = min(len(s1),len(s2))
  shared = 0
  for ind in range(minLength-1):
    if s1[ind] == s2[ind]: shared+=1
    else: break
  if minLength == 1:
    return max(0,shared-1)
  else:
    return shared

def sharedSuffix(s1,s2):
  minLength = min(len(s1),len(s2))-1
  shared = 0
  for ind in range(minLength*-1,0)[::-1]:
    if s1[ind] == s2[ind]: shared+=1
    else: break
  return shared

###################
# VCF FIELDS
###################

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

outsnvs = sys.stdout
outindels = sys.stdout
if options.SNVs != None:
  outsnvs = open(options.SNVs,'w')
if options.InDels != None:
  outindels = open(options.InDels,'w')

for line in sys.stdin:
  line = line.rstrip('\r\n')
  if line.upper().startswith("#CHROM"):
    fields = line.split()
    outsnvs.write("\t".join(fields[:falt_allele+1])+"\n")
    if options.SNVs != None or options.InDels != None:
      outindels.write("\t".join(fields[:falt_allele+1])+"\n")
  if line.startswith("#"): continue
  fields = line.split("\t")
  if len(fields) > falt_allele and is_nucleotide(fields[fref_allele]):
    ref = fields[fref_allele].upper()
    if options.refCheck:
      pos = int(fields[fpos])
      href = reference.fetch(fields[fchr],pos-1,pos-1+len(ref)).upper()
      if href != ref: continue
      
    for alt in fields[falt_allele].upper().split(','):
      if ref == alt: continue
      if is_nucleotide(alt):
        if len(alt) == len(ref) and len(ref) == 1:
          outsnvs.write("%s\t%s\t.\t%s\t%s\n"%(fields[fchr],fields[fpos],ref,alt))
        else:
          trimValue = sharedPrefix(ref,alt)
          if trimValue != 0:
            nref = ref[trimValue:]
            nalt = alt[trimValue:]
          else:
            nref,nalt = ref,alt
          trimValue2 = sharedSuffix(nref,nalt)
          if trimValue2 != 0:
            nref = nref[:-trimValue2]
            nalt = nalt[:-trimValue2]
          if len(nalt) == len(nref) and len(ref) == 1:
            outsnvs.write("%s\t%d\t.\t%s\t%s\n"%(fields[fchr],int(fields[fpos])+trimValue,nref,nalt))
          elif (trimValue == 0):
            outindels.write("%s\t%s\t.\t%s\t%s\n"%(fields[fchr],fields[fpos],nref,nalt))
          else:
            outindels.write("%s\t%d\t.\t%s\t%s\n"%(fields[fchr],int(fields[fpos])+trimValue,nref,nalt))

if options.SNVs != None:
  outsnvs.close()
if options.InDels != None:
  outindels.close()
