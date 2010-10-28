#!/usr/bin/python
"""
  Date of Creation: 3rd Apr 2010    
                       
  Description:   Classes and functions for manipulating BED files

  Copyright (C) 2010  
  University of Southern California,
  Philip J. Uren,
  Jin H. Park,
  Andrew D. Smith
  
  Authors: Philip J. Uren, Jin H. Park, Andrew D. Smith
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
  --------------------
  
  Known Bugs:    None
  
  Revision 
  History:       18th August 2010 -- Philip Uren
                   * added BEDIterator 
                   * updated header details 
                   * added testSizeOfOverlap unit test
                 20th August 2010 -- Philip Uren
                   * added check for same chrom to intersection test
                 13th September 2010 -- Philip Uren
                   * modified BEDElement and BEDElementFromString to allow 
                     BED elements to have only chrom, start and end 
                   * cleaned details of junction reads
                 16th September 2010 -- Philip Uren
                   * added check for sorted order to BEDIterator 
                 17th September 2010 -- Philip Uren
                   * added intervalTree function 
                 24th September 2010 -- Philip Uren
                   * fixed bugs with checking for sorted BED files when 
                     not required 
                 1st October 2010 -- Philip Uren
                   * added ability to check BED file is sorted by read name
                   * added BEDUniqueFilter
                   * added unit tests for above two items
                   * added missing exception class 
                   * added unit test for BED file being sorted by start
                   * added verbose option to BEDIterator
                 2nd October 2010 -- Philip Uren
                   * added BEDDuplicateFilter and associated unit tests
                 19th October 2010 -- Philip Uren
                   * added best option for BEDUniqueFilter
                   * changed type for score from float to int 
                 22nd October 2010 -- Philip Uren
                   * fixed bug in BEDIterator when passing string instead of
                     stream 
                   * added detail to the exception raised when BED elements 
                     appear not to have enough parts to them 
                 23rd October 2010 -- Philip Uren
                   * added BEDIterator sorting order of chromosome 
                   * added unit test for above 
                 27th October 2010 -- Philip Uren
                   * moved to smithlab_py 
                   * changed usage of progress indicator to use new class
                 28th October 2010 -- Philip Uren 
                   * changed unit test for overlap
                   * removed iterators to a new module 
  
  TODO:         
                * finish unit tests

"""

import sys
import os
import unittest
import copy
from Queue import Queue

from testing.dummyfiles import DummyInputStream, DummyOutputStream
from util.progressIndicator import ProgressIndicator
from util.fileUtils import linesInFile
from datastruct.intervalTree import IntervalTree


DEFAULT_DELIM = "\t"
ITERATOR_SORTED_START = 1
ITERATOR_SORTED_END = 2
ITERATOR_SORTED_NAME = 3
ITERATOR_SORTED_CHROM = 4

class BEDError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)

def intervalTrees(reffh, verbose = False):
  """
    DESCRP: Build a dictionary of interval trees indexed by chrom from
            a BED stream
    PARAMS: reffh -- a stream allowing iteration over lines in BED format,
                     or a filename for a BED file
            verbose -- output progress messages to sys.stderr if True
  """
  if type(reffh).__name__ == "str" : fh = open(reffh)
  else : fh = reffh
  
  # load all the regions and split them into lists for each chrom
  elements = {}
  if verbose and fh != sys.stdin :
    totalLines = linesInFile(fh.name)
    pind = ProgressIndicator(totalToDo = totalLines, 
                                   messagePrefix = "completed", 
                                   messageSuffix = "of loading " + fh.name)      
  for element in BEDIterator(fh):
    if not element.chrom in elements : elements[element.chrom] = []
    elements[element.chrom].append(element)
    if verbose and fh != sys.stdin:
      pind.done += 1
      pind.showProgress()
    
  # create an interval tree for each list
  trees = {}
  if verbose:
    totalLines = len(elements)
    pind = ProgressIndicator(totalToDo = totalLines, 
                                   messagePrefix = "completed", 
                                   messageSuffix = "of making interval trees")
  for chrom in elements :
    trees[chrom] = IntervalTree(elements[chrom])
    if verbose:
      pind.done += 1
      pind.showProgress()
      
  return trees


def BEDElementFromString(line, delim = DEFAULT_DELIM):
  peices = line.split(delim)
  if len(peices) < 3 : 
    raise BEDError("BED elements must have at least chrom, start and end" +\
                   " found only " + str(len(peices)) + " in " + line)
  chrom = peices[0]
  start = peices[1]
  end = peices[2]
  
  name = None 
  score = None 
  strand = None 
  
  if len(peices) >= 4 != None : name = peices[3]
  if len(peices) >= 5 != None : score = peices[4]
  if len(peices) >= 6 != None : strand = peices[5]
  
  return BEDElement(chrom, start, end, name, score, strand)
  

class BEDElement :
  def __init__(self, chrom, start, end, name, score, strand, delim = DEFAULT_DELIM):
    """
      ...
    """
    # the basic read info
    self.chrom = chrom.strip()
    self.start = int(start)
    self.end = int(end)
    
    # we might get the following too...
    self.name = None
    self.score = None
    self.strand = None    
    if name != None: self.name = name.strip()
    if score != None : self.score = int(score)
    if strand != None : self.strand = strand.strip()
    
    # we use this to keep track of the reads that map to this exon
    self.readsMapped = []
    
    # use this to keep track of whether this read has mapped to any genes yet
    self.mapped = False
    
    # just remember what delimiters we used, so we can re-construct the string later
    self.delim = delim 
  
  def __eq__(self, e):
    if e == None : return False
    return self.chrom == e.chrom and self.start == e.start and\
           self.end == e.end and self.name == e.name and\
           self.score == e.score and self.strand == e.strand
        
  def __len__(self):
    return (self.end - self.start) + 1      
        
  def __str__(self):
    res = self.chrom + self.delim + str(self.start) + self.delim + str(self.end)
    if self.name != None : res += (self.delim + self.name)
    if self.score != None : res += (self.delim + str(self.score))
    if self.strand != None : res += (self.delim + self.strand)
     
    return  res
                
  def readCountsAsString(self):
    return self.name + " = " +str(self.readsMapped)
  
  def distance(self, e):
    dist = 0
    if not e.intersects(self) :
      dist = max(self.start, e.start) - min(self.end, e.end)
    return dist
  
  def subtract(self, e):
    """ 
      subtracts the region defined by e from this region -- if e is entirely inside self, then e is treated as extending to 
      the closest of either start or end. If self is entirely inside e, then we arbitrarily set self to size 0 by setting end = start 
    """
    es = e.start
    ee = e.end
    
    # check that they overlap...
    if not self.intersects(e) :
      return
    
    # handle case where e is inside self
    if es >= self.start and es <= self.start and ee >= self.start and ee <= self.end :
      sdiff = math.abs(self.start - es)
      ediff = math.abs(self.end - ee)
      if sdiff > ediff : ee = self.end
      else : es = self.start
    
    # handle case where self is inside e 
    if self.start <= ee and self.start >= es and self.end >= es and self.end <= ee :
      self.end = self.start
      return
    
    # general case - some overlap by neither includes the other entirely 
    toShrinkEnd = self.end - es 
    toShrinkStart = ee - self.start
    if toShrinkEnd < toShrinkStart : self.end = es
    else : self.start = ee
    
  def sizeOfOverlap(self, e):
    # no overlap
    if not self.intersects(e) : return 0
    
    # complete inclusion..
    if e.start >= self.start and e.end <= self.end : return len(e)
    if self.start >= e.start and self.end <= e.end : return len(self)
    
    # partial overlap
    if e.start > self.start : return (self.end - e.start + 1)
    if self.start > e.start : return (e.end - self.start + 1)
    
    
  def intersects(self, e):
    """ returns true if this elements intersects the element e """

    if self.chrom != e.chrom : return False
    if self.end >= e.start and self.end <= e.end : return True
    if self.start >= e.start and self.start <= e.end : return True
    if e.start >= self.start and e.start <= self.end : return True
    if e.end >= self.start and e.end <= self.end : return True
    return False
    


class BEDUnitTests(unittest.TestCase):
  """
    Unit tests for fqsplit 
  """
  
  def setUp(self):
    pass
    
  def testInterection(self):
    pass
    
      
  def testSizeOfOverlap(self):
      ## region 2 entirely inside region 1 -- ans = size of region 2
      region1 = BEDElement("chr1", 10, 20, "read", 1, "f")
      region2 = BEDElement("chr1", 10, 20, "read", 1, "f")
      self.assertEqual(region1.sizeOfOverlap(region2), len(region2))
      self.assertEqual(region1.sizeOfOverlap(region2), 
                       region2.sizeOfOverlap(region1))
  
  def testSubtract(self):
    pass
  
  def testDistance(self):
    pass
  

if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])
