#!/usr/bin/env python
# -*- coding: ASCII -*-

import sys, os
import pysam
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p","--path", dest="path", help="Path to scored variants.")
parser.add_option("-i","--input", dest="input", help="Read variants from vcf file (default stdin)",default=None)
parser.add_option("--found_out", dest="found_out", help="Write found variants to file (default: stdout)",default=None)
parser.add_option("--header", dest="header", help="Write full header to output (default none)",default=False, action="store_true")
(options, args) = parser.parse_args()

if options.input:
    stdin = open(options.input,'r')
else:
    stdin = sys.stdin

if options.found_out:
    found_out = open(options.found_out,'w')
else:
    found_out = sys.stdout

fpos,fref,falt = 1,2,3
if os.path.exists(options.path) and os.path.exists(options.path+".tbi"):
  filename = options.path
  sys.stderr.write("Opening %s...\n"%(filename))
  regionTabix = pysam.Tabixfile(filename,'r')
  header = list(regionTabix.header)
  for line in header:
    if options.header:
      found_out.write(line+"\n")
    try:
      fref = line.split('\t').index('Ref')
      falt = line.split('\t').index('Alt')
    except ValueError:
      pass
else:
  raise IOError("No valid file with pre-scored variants.\n")

for line in stdin:
  line = line.rstrip('\n\r')
  if line.startswith('#'):
    sys.stdout.write(line + '\n')
    continue

  try:
    fields = line.split('\t')
    found = False
    chrom = fields[0]
    pos = int(fields[1])
    lref,allele = fields[-2],fields[-1]
    for regionHit in regionTabix.fetch(chrom, pos-1, pos):
      vfields = regionHit.rstrip().split('\t')
      if (vfields[fref] == lref) and (vfields[falt] == allele) and (vfields[fpos] == pos)):
        found_out.write(regionHit+"\n")
        found = True
        break

    if not found:
      sys.stdout.write(line + '\n')

  except ValueError:
    sys.stderr.write('Encountered uncovered chromosome\n')
    sys.stdout.write(line + '\n')

if options.input:
    stdin.close()

if options.found_out:
    found_out.close()
