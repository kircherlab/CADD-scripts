#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *27.03.2008

"""

import math

class EDict:

  def __init__(self):
    self._internal_dict = {}
    self._sorted_keys = []
    self._is_sync = True

  def __len__(self):
    return len(self._internal_dict.keys())

  def __getitem__(self,key):
    return self._internal_dict[key]

  def __repr__(self):
      """Return string representation of a EDict."""
      return repr(self._internal_dict)

  # __str__ is the same as __repr__
  __str__ = __repr__

  def clear(self):
    self._is_sync = True
    self._internal_dict = {}
    self._sorted_keys = []

  def set(self,key,value):
    if key == None:
      raise TypeError,"EDict does not allow None keys"
    if not(self._internal_dict.has_key(key)):
      self._is_sync = False
    self._internal_dict[key] = value
    return True


  def get(self,key):
    if not(self._internal_dict.has_key(key)):
      return None
    else:
      return self._internal_dict[key]

  def get_keys(self): return self._internal_dict.keys()

  def get_values(self): return self._internal_dict.values()

  def get_items(self): return self._internal_dict.items()

  def iterkeys(self): return self._internal_dict.iterkeys()

  def itervalues(self): return self._internal_dict.itervalues()

  def iteritems(self): return self._internal_dict.iteritems()

  def has_key(self,elem): return self._internal_dict.has_key(elem)

  def _make_sync(self):
    self._sorted_keys = self._internal_dict.keys()
    self._sorted_keys.sort()
    self._is_sync = True

  def get_smaller(self,key):
    if not(self._is_sync): self._make_sync()
    cur_len = len(self._sorted_keys)
    if cur_len > 0:
      if not(key <= self._sorted_keys[0]):
        cur_pos = -1
        forlast = -1
        new_pos = cur_len/2
        dist = max(int(round(cur_len/4.0)),1)
        while(cur_pos != new_pos) and (new_pos != forlast):
          forlast = cur_pos
          cur_pos = new_pos
          #print cur_pos, dist
          if key > self._sorted_keys[cur_pos]:
            new_pos = cur_pos+dist
            if (new_pos >= cur_len): new_pos = cur_len-1
          elif key < self._sorted_keys[cur_pos]:
            new_pos = cur_pos-dist
            if (new_pos < 0): new_pos = 0
          else:
            new_pos = cur_pos
          dist = max(int(dist/2.0),1)
        if (cur_pos+1 < cur_len) and (cur_pos > 0):
          if (key == self._sorted_keys[cur_pos]):
            return (self._sorted_keys[cur_pos-1],self._internal_dict[self._sorted_keys[cur_pos-1]])
          elif ((key <= self._sorted_keys[cur_pos+1]) and (key > self._sorted_keys[cur_pos])):
            return (self._sorted_keys[cur_pos],self._internal_dict[self._sorted_keys[cur_pos]])
          elif ((key <= self._sorted_keys[cur_pos]) and (key > self._sorted_keys[cur_pos-1])):
            return (self._sorted_keys[cur_pos-1],self._internal_dict[self._sorted_keys[cur_pos-1]])
          else:
            print "get_smaller: SHOULD NOT HAPPEN!",cur_pos,"max:",cur_len
        elif ((cur_pos == 0) and (key > self._sorted_keys[cur_pos])):
          return (self._sorted_keys[cur_pos],self._internal_dict[self._sorted_keys[cur_pos]])
        elif ((cur_pos+1 == cur_len) and (key > self._sorted_keys[cur_pos])):
          return (self._sorted_keys[cur_pos],self._internal_dict[self._sorted_keys[cur_pos]])
        elif ((cur_pos+1 == cur_len) and (key > self._sorted_keys[cur_pos-1])):
          return (self._sorted_keys[cur_pos-1],self._internal_dict[self._sorted_keys[cur_pos-1]])
    return (None,None)

  def get_smaller_equal(self,key):
    if self.has_key(key):
      return key,self._internal_dict[key]
    else:
      return self.get_smaller(key)

  def get_larger(self,key):
    if not(self._is_sync): self._make_sync()
    cur_len = len(self._sorted_keys)
    if cur_len > 0:
      if not(key >= self._sorted_keys[-1]):
        cur_pos = -1
        forlast = -1
        new_pos = cur_len/2
        dist = max(int(round(cur_len/4.0)),1)
        while(cur_pos != new_pos) and (new_pos != forlast):
          forlast = cur_pos
          cur_pos = new_pos
          #print cur_pos, dist
          if key > self._sorted_keys[cur_pos]:
            new_pos = cur_pos+dist
            if (new_pos >= cur_len): new_pos = cur_len-1
          elif key < self._sorted_keys[cur_pos]:
            new_pos = cur_pos-dist
            if (new_pos < 0): new_pos = 0
          else:
            new_pos = cur_pos
          dist = max(int(dist/2.0),1)
        if (cur_pos+1 < cur_len) and (cur_pos > 0):
          if (key == self._sorted_keys[cur_pos]):
            return (self._sorted_keys[cur_pos+1],self._internal_dict[self._sorted_keys[cur_pos+1]])
          elif ((key < self._sorted_keys[cur_pos+1]) and (key >= self._sorted_keys[cur_pos])):
            return (self._sorted_keys[cur_pos+1],self._internal_dict[self._sorted_keys[cur_pos+1]])
          elif ((key < self._sorted_keys[cur_pos]) and (key >= self._sorted_keys[cur_pos-1])):
            return (self._sorted_keys[cur_pos],self._internal_dict[self._sorted_keys[cur_pos]])
          else:
            print "get_larger: SHOULD NOT HAPPEN!",cur_pos,"max:",cur_len
        elif ((cur_pos == 0) and (key < self._sorted_keys[cur_pos])):
          return (self._sorted_keys[cur_pos],self._internal_dict[self._sorted_keys[cur_pos]])
        elif ((cur_pos == 0) and (key < self._sorted_keys[cur_pos+1])):
          return (self._sorted_keys[cur_pos+1],self._internal_dict[self._sorted_keys[cur_pos+1]])
        elif ((cur_pos+1 == cur_len) and (key < self._sorted_keys[cur_pos])):
          return (self._sorted_keys[cur_pos],self._internal_dict[self._sorted_keys[cur_pos]])
    return (None,None)

  def get_larger_equal(self,key):
    if self.has_key(key):
      return key,self._internal_dict[key]
    else:
      return self.get_larger(key)


if __name__ == '__main__':
  import os
  
  convtbl = EDict()
  convtbl.set(149.449379,"99")
  convtbl.set(149.399739,"96")
  convtbl.set(4.393534,"23")
  convtbl.set(3,"15")
  convtbl.set(0,"5")
  convtbl.set(-1.933954,"0.007")
  convtbl.set(-2.118821,"0.005")
  convtbl.set(-2.388300,"0.003")
  convtbl.set(-2.874310,"0.001")
  convtbl.set(-226.479313,"0")
  print 150,convtbl.get_larger_equal(150)
  print 149.4,convtbl.get_larger_equal(149.4)
  print 3,convtbl.get_larger_equal(3)
  print -2.3,convtbl.get_larger_equal(-2.3)
  print -2.4,convtbl.get_larger_equal(-2.4)
  print -2.8,convtbl.get_larger_equal(-2.8)
  print -2.9,convtbl.get_larger_equal(-2.9)
  print -227,convtbl.get_larger_equal(-227)

  print "###############"
  convtbl = EDict()
  convtbl.set(149.449379,"99")
  convtbl.set(149.399739,"96")
  convtbl.set(4.393534,"23")
  convtbl.set(3,"15")
  convtbl.set(-1.933954,"0.007")
  convtbl.set(-2.018751,"0.006")
  convtbl.set(-2.118821,"0.005")
  convtbl.set(-2.239310,"0.004")
  convtbl.set(-2.388300,"0.003")
  convtbl.set(-2.580755,"0.002")
  convtbl.set(-2.874310,"0.001")
  convtbl.set(-226.479313,"0")
  print 150,convtbl.get_larger_equal(150)
  print 149.4,convtbl.get_larger_equal(149.4)
  print 3,convtbl.get_larger_equal(3)
  print -2.3,convtbl.get_larger_equal(-2.3)
  print -2.4,convtbl.get_larger_equal(-2.4)
  print -2.6,convtbl.get_larger_equal(-2.6)
  print -2.8,convtbl.get_larger_equal(-2.8)
  print -2.9,convtbl.get_larger_equal(-2.9)
  print -227,convtbl.get_larger_equal(-227)
  
  print "###############"
  table = os.environ['CADD'] + "/whole_genome/conversion_tbl_cave/conversion_table_ext.tsv"
  maxValue,minValue = None,None
  convtbl = EDict()
  if os.path.exists(table):
    infile = open(table)
    for line in infile:
      fields = line.split()
      if len(fields) == 2:
        val = float(fields[1])
        convtbl.set(val,fields[0])
        if val > maxValue or maxValue == None: maxValue = val
        if val < minValue or minValue == None: minValue = val
    infile.close()
    #convtbl.set(-220.0,"0")
    print 150,convtbl.get_larger_equal(150)
    print 149.4,convtbl.get_larger_equal(149.4)
    print 3,convtbl.get_larger_equal(3)
    print -2.3,convtbl.get_larger_equal(-2.3)
    print -2.4,convtbl.get_larger_equal(-2.4)
    print -2.6,convtbl.get_larger_equal(-2.6)
    print -2.8,convtbl.get_larger_equal(-2.8)
    print -2.9,convtbl.get_larger_equal(-2.9)
    print -227,convtbl.get_larger_equal(-227)
    
    print "###"
    print len(convtbl)
    print "###"
    count = 0
    for key,value in sorted(convtbl.iteritems()):
      print key,value
      count += 1
      if count > 10: break
