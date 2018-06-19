#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *02.12.2012

"""


import sys, os
from lib import EDict
import math
from optparse import OptionParser

parser = OptionParser("%prog [options]")
parser.add_option("-t", "--table", dest="table", help="Scoring table (default conversion_table_ext.tsv)", default="conversion_table_ext.tsv")
parser.add_option("-r", "--replace", dest="replace", help="Files are already scored, replace the score value (def OFF)", default=False, action="store_true")
(options, args) = parser.parse_args()

maxValue,minValue = None,None
convtbl = EDict.EDict()
if os.path.exists(options.table):
  infile = open(options.table)
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      val = float(fields[1])
      convtbl.set(val,fields[0])
      if val > maxValue or maxValue == None: maxValue = val
      if val < minValue or minValue == None: minValue = val
  infile.close()
 
replace = options.replace
for line in sys.stdin:
  if line.startswith('#') or line.startswith("Chrom"):
    if "PHRED" not in line:
      sys.stdout.write(line.rstrip()+"\tPHRED\n")
    else:
      sys.stdout.write(line)
      replace = True
  else:
    fields = line.rstrip().split('\t')
    if replace:
      val = float(fields[-2])
      fields[-2]="%.6f"%(val)
    else:
      val = float(fields[-1])
      fields[-1]="%.6f"%(val)
    score = convtbl.get_larger_equal(val)
    if score[-1] == None:
      if val > maxValue:
        score = convtbl.get_larger_equal(maxValue)
      elif val < minValue:
        score = convtbl.get_larger_equal(minValue)
      else:
        score = (None,"NA")
        
    if replace:
      fields[-1] = score[-1]
      sys.stdout.write("\t".join(fields)+"\n")
    else:
      sys.stdout.write("\t".join(fields)+"\t%s\n"%(score[-1]))
