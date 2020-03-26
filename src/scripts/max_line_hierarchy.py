#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *06.07.2013
"""

import sys, os
from optparse import OptionParser
import gzip

parser = OptionParser("%prog [options]")
parser.add_option("--noheader", dest="noheader", help="Do not print header",default=False,action="store_true")
parser.add_option("--PHRED", dest="PHRED", help="File is already PHRED-scaled",default=False,action="store_true")
parser.add_option("-a","--all", dest="all", help="Report all annotation lines of variant with max score",default=False,action="store_true")
parser.add_option("-i","--input", dest="input", default=None, type=str, help="Input tsv.gz file (default=stdin)")
(options, args) = parser.parse_args()

if options.input:
    stdin = gzip.open(options.input)
else:
    stdin = sys.stdin

if options.PHRED:
  scorecol = -2
else:
  scorecol = -1

fVconseq = 7
annolines = [(None,None),(None,None),(None,None),(None,None)]
all_lines = []

lpos,maxline,maxscore = None,None,None

# retain GPS values to the highest scoring variant if the one having them is not the one being kept
GPSvalues = ["NA"]
GPSindices = range(109, 115)
GPSnames = ['Grantham', 'PolyPhenCat', 'PolyPhenVal','SIFTcat','SIFTval']

for line in stdin:
  if line.startswith('#'):
    fields = line.rstrip().split("\t")
    GPSindices = [ind for ind,value in enumerate(fields) if value in GPSnames]
    for ind, value in enumerate(fields):
        if value == 'RawScore':
            scorecol = ind
            break
        if value == 'Consequence':
            fVconseq = ind
    if options.noheader: continue
    sys.stdout.write(line)
  else:
    fields = line.rstrip().split('\t')
    score = float(fields[scorecol])
    pos = tuple(fields[:5])

    if lpos != pos:
      if lpos != None:
        if not options.all:
          for maxscore,maxline in annolines:
            if maxscore != None:
              if (maxline[GPSindices[0]] == "NA") and (GPSvalues[0] != "NA"):
                for ind, val in zip(GPSindices, GPSvalues):
                    maxline[ind] = val
              sys.stdout.write("\t".join(maxline)+"\n")
              break
        else:
          replace = None
          for maxscore,maxline in annolines:
            if maxscore != None:
              replace = maxline[scorecol:]
              break
          for linefields in all_lines:
            sys.stdout.write("\t".join(linefields[:scorecol]+replace)+"\n")
          all_lines = []
        
      GPSvalues = [fields[ind] for ind in GPSindices]
      annolines = [(None,None),(None,None),(None,None),(None,None)]
      lpos = pos

    if fields[fVconseq] in ['NON_SYNONYMOUS', 'SYNONYMOUS', 'STOP_GAINED', 'STOP_LOST', 'FRAME_SHIFT', 'INFRAME', 'CANONICAL_SPLICE', 'SPLICE_SITE', 'NONCODING_CHANGE']:
      if annolines[0][0] == None or score > annolines[0][0]:
        annolines[0]=score,fields
    elif fields[fVconseq] in ['REGULATORY', '5PRIME_UTR', '3PRIME_UTR']:
      if annolines[1][0] == None or score > annolines[1][0]:
        annolines[1]=score,fields
    elif fields[fVconseq] in ['INTRONIC', 'UPSTREAM', 'DOWNSTREAM']:
      if annolines[2][0] == None or score > annolines[2][0]:
        annolines[2]=score,fields
    else:
      if annolines[3][0] == None or score > annolines[3][0]:
        annolines[3]=score,fields

    if (GPSvalues[0] == "NA") and (fields[GPSindices[0]] != "NA"):
      GPSvalues = [fields[ind] for ind in GPSindices]

    if options.all:
      all_lines.append(fields)

if lpos != None:
  if not options.all:
    for maxscore,maxline in annolines:
      if maxscore != None:
        if (maxline[GPSindices[0]] == "NA") and (GPSvalues[0] != "NA"):
          for ind, val in zip(GPSindices, GPSvalues):
            maxline[ind] = val
        sys.stdout.write("\t".join(maxline)+"\n")
        break
  else:
    replace = None
    for maxscore,maxline in annolines:
      if maxscore != None: 
        replace = maxline[scorecol:]
        break
    for linefields in all_lines:
      sys.stdout.write("\t".join(linefields[:scorecol]+replace)+"\n")
    all_lines = []
