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

def BEDDuplicateFilter(fh1, fh2, removeJuncTags=False, removePETags=False, 
                       verbose=False):
  """
    @summary: returns an iterator which yields pairs of reads from fh1
              and fh2 such that duplicate reads (same name) are included
              in the same pair, while unique reads are paired with None.
              Reads from fh1 will always be first in the tuple, while reads
              from fh2 will always be second.
    @param fh1: BED formated stream, must be sorted by name
    @param fh2: BED formated stream, must be sorted by name
    @param removeJuncTags: if True, remove the junctions tags added to the 
                           read names by juncConvert (if present) before 
                           comparing names (output reads will still have them)
    @param removePETags: if True, remove the paired end tags of the read names
                         before comparing. 
    @param verbose: output additional progress messages to stderr if this 
                    is True         
    @note: it is assumed that PE tags are .1 or .2 and junc tags are %1 or %2
           and that junc tags will appear after PE tags
  """
  
  if removePETags and not removeJuncTags :
    raise BEDError("BEDDuplicateFilter cannot trim PE tags if junc tags are" +\
                   "not also trimmed -- maybe in a future version")
  
  def next(iterator):
    """ little internal function to return the next item, or None """
    try : return iterator.next()
    except StopIteration : return None
  
  def name(read):
    name = read.name
    if removeJuncTags and name[-2] == "%" : name = name[:-2]  
    if removePETags and name[-2] == "." : name = name[:-2]
    return name
  
  def match(read, list1, list2):
    """ 
      return True if the name in read can be added to the lists without
      destroying the homogeneity of the names in the lists  
    """
    for r in list1 : 
      if name(r) != name(read) : return False 
    for r in list2 : 
      if name(r) != name(read) : return False
    return True
  
  bit1 = BEDIterator(fh1, sortedby=ITERATOR_SORTED_NAME, verbose=verbose)
  bit2 = BEDIterator(fh2, sortedby=ITERATOR_SORTED_NAME, verbose=verbose)
  bf1, bf2 = None, None
  bf1List, bf2List = [], []
  
  while True:
    # try to get an item from each stream if we don't already have one
    if bf1 == None: bf1 = next(bit1)
    if bf2 == None: bf2 = next(bit2)
    
    if bf1 != None and (bf2 == None or name(bf1) < name(bf2)) :
      if not match(bf1, bf1List, bf2List) :
        yield (bf1List, bf2List)
        bf1List, bf2List = [], []
      bf1List.append(bf1)
      bf1 = next(bit1)
    elif bf2 != None and (bf1 == None or name(bf2) < name(bf1)) :
      if not match(bf2, bf1List, bf2List) :
        yield (bf1List, bf2List)
        bf1List, bf2List = [], []
      bf2List.append(bf2)
      bf2 = next(bit2)
    elif bf1 != None and bf2 != None and name(bf1) == name(bf2) :
      if not match(bf2, bf1List, bf2List) :
        yield (bf1List, bf2List)
        bf1List, bf2List = [], []
      bf1List.append(bf1)
      bf2List.append(bf2)
      bf1 = next(bit1)
      bf2 = next(bit2)
      
    if bf1 == None and bf2 == None:
      if bf1List != [] or bf2List != None :
        yield (bf1List, bf2List)
        bf1List, bf2List = [], []
      break

def BEDUniqueFilter(fh1, fh2, verbose=False, best=False):
  """
    @summary: returns an iterator which will yield pairs of reads (r1, r2)
              for which r1 does not appear in fh2 and r2 does not appear in
              fh1. Either r1 or r2 may be None (both not both) if one of 
              file1 or file2 has fewer unique entries than the other. 
              reads are compared by name, other attributes are ignored 
    @param fh1: BED formated stream, must be sorted by name
    @param fh2: BED formated stream, must be sorted by name
    @param verbose: output additional progress messages to stderr if this 
                    is True 
    @param best: If True, when encountering two reads with the same name, 
                 we'll output the one with the better (smaller) score, unless
                 they both get the same score, then we'll skip them
  """
  
  def next(iterator):
    """ little internal function to return the next item, or None """
    try : return iterator.next()
    except StopIteration : return None   
  
  bit1 = BEDIterator(fh1, sortedby=ITERATOR_SORTED_NAME, verbose=verbose)
  bit2 = BEDIterator(fh2, sortedby=ITERATOR_SORTED_NAME, verbose=verbose)
  bf1, bf2 = None, None
  bf1_Q, bf2_Q = Queue(), Queue()
  bf1_exhausted, bf2_exhausted = False, False
  
  while True:
    # try to get an item from each stream if we don't already have one
    if bf1 == None: bf1 = next(bit1)
    if bf2 == None: bf2 = next(bit2)
    
    # the item in bf1 (if there is one) can be output if it's less
    # than bf2, or we've used up all of the second stream
    if bf1 != None and (bf2 == None or bf1.name < bf2.name) :
      bf1_Q.put(bf1)
      bf1 = next(bit1)
    if bf2 != None and (bf1 == None or bf2.name < bf1.name) :
      bf2_Q.put(bf2)
      bf2 = next(bit2)
    if bf1 != None and bf2 != None and bf1.name == bf2.name :
      if best and bf1.score < bf2.score : bf1_Q.put(bf1)
      if best and bf2.score < bf1.score : bf2_Q.put(bf2)
      bf1 = next(bit1)
      bf2 = next(bit2)
      
    # output a pair
    y1, y2 = None, None
    if not bf1_Q.empty() : y1 = bf1_Q.get()
    if not bf2_Q.empty() : y2 = bf2_Q.get()
    yield (y1,y2)
    
    if bf1_Q.empty() and bf2_Q.empty() and bf1 == None and bf2 == None:
      break
    
    

def BEDIterator(filehandle, sortedby=None, verbose=False):
  """
    DESCRP: Get an iterator for a BED file
    PARAMS: sortedby -- if None, order is not checked.
                        if == ITERATOR_SORTED_START, elements in file must
                          be sorted by chrom and start index (an exception 
                          is raised if they are not)
                        if == ITERATOR_SORTED_END, element must be sorted
                          by chrom and end index
    RETURN: iterator where subsequent calls to next() yield the next element
            in the file
  """
  chromsSeen = set()
  prev = None
  if type(filehandle).__name__ == "str" : filehandle = open(filehandle)
  
  if verbose :
    try :
      totalLines = linesInFile(filehandle.name)
      pind = ProgressIndicator(totalToDo = totalLines, 
                                     messagePrefix = "completed", 
                                     messageSuffix = "of processing " +\
                                                      filehandle.name)
    except AttributeError :
      sys.stderr.write("BEDIterator -- warning: " +\
                       "unable to show progress for stream")
      verbose = False      

  
  for line in filehandle :
    if verbose :
      pind.done += 1
      pind.showProgress()
      
    if line.strip() == "": continue
    e = BEDElementFromString(line)
    
    # sorting by name?
    if sortedby == ITERATOR_SORTED_NAME and prev != None and prev.name > e.name :
      raise BEDError("bed file " + filehandle.name +\
                       " not sorted by element name" +\
                       " found " + e.name + " after " +\
                       prev.name)
    
    # first item
    if prev == None :
      chromsSeen.add(e.chrom)
    
    # on same chrom as the prev item, make sure order is right
    if prev != None and sortedby != None and e.chrom == prev.chrom :
      if sortedby == ITERATOR_SORTED_START and prev.start > e.start :
        raise BEDError("bed file " + filehandle.name +\
                       " not sorted by start index - saw item " +\
                       str(prev) + " before " + str(e))
      if sortedby == ITERATOR_SORTED_END and prev.end > e.end :
        raise BEDError("bed file " + filehandle.name +\
                       " not sorted by end index")
    
    # starting a new chrom.. make sure we haven't already seen it
    if prev != None and prev.chrom != e.chrom :
      if (sortedby == ITERATOR_SORTED_START or 
          sortedby == ITERATOR_SORTED_END or
          sortedby == ITERATOR_SORTED_CHROM) and\
         (e.chrom in chromsSeen or prev.chrom > e.chrom) :
        raise BEDError("BED file " + filehandle.name +\
                       " not sorted by chrom")
      chromsSeen.add(e.chrom) 
    
    # all good..
    yield e
    prev = e

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

  def testBEDUniqueFilter(self):
    debug = False
    end1 = "chr1 \t 10 \t 14 \t read1\n" +\
                "chr1 \t 30 \t 34 \t read2\n" +\
                "chr1 \t 11 \t 15 \t read3"
    end2 = "chr1 \t 15 \t 19 \t read1\n" +\
                "chr1 \t 80 \t 84 \t read2\n" +\
                "chr1 \t 80 \t 84 \t read4\n" +\
                "chr1 \t 80 \t 84 \t read5\n"
    keep1 = ["chr1\t11\t15\tread3"]
    keep2 = ["chr1\t80\t84\tread4",
             "chr1\t80\t84\tread5"]
    keep1.sort()
    keep2.sort()
    
    end1_infs = DummyInputStream(end1)
    end2_infs = DummyInputStream(end2)
    gotOutput1 = []
    gotOutput2 = []
    for r1,r2 in BEDUniqueFilter(end1_infs, end2_infs):
      if r1 != None : gotOutput1.append(str(r1))
      if r2 != None : gotOutput2.append(str(r2))
                                        
    gotOutput1.sort()
    gotOutput2.sort()
    
    if debug :
      print "got for 1 ------ "
      for l in gotOutput1 : print l
      print "-------"
      print "got for 2 ------ "
      for l in gotOutput2 : print l
      print "-------"
    
    self.assertTrue(gotOutput1 == keep1)
    self.assertTrue(gotOutput2 == keep2)
    
  def testBEDDuplicateFilter(self):
    debug = False
    end1 = "chr1 \t 10 \t 14 \t read1.1\n" +\
           "chr1 \t 30 \t 34 \t read3.1\n" +\
           "chr1 \t 11 \t 15 \t read5.1\n" +\
           "chr1 \t 11 \t 15 \t read8.1\n" +\
           "chr1 \t 11 \t 15 \t read9.1%1\n" +\
           "chr1 \t 11 \t 15 \t read9.1%2\n"
    end2 = "chr1 \t 15 \t 19 \t read2.2\n" +\
           "chr1 \t 80 \t 84 \t read3.2%1\n" +\
           "chr1 \t 80 \t 84 \t read3.2%2\n" +\
           "chr1 \t 80 \t 84 \t read3.2%3\n" +\
           "chr1 \t 80 \t 84 \t read4.2\n" +\
           "chr1 \t 80 \t 84 \t read5.2\n" +\
           "chr1 \t 80 \t 84 \t read9.2\n"
    expectOutput = [(["read1.1"], []),
                   ([], ["read2.2"]),
                   (["read3.1"], ["read3.2%1","read3.2%2", "read3.2%3"]),
                   ([], ["read4.2"]),
                   (["read5.1"], ["read5.2"]),
                   (["read8.1"],[]),
                   (["read9.1%1","read9.1%2"],["read9.2"])]
    
    end1_infs = DummyInputStream(end1)
    end2_infs = DummyInputStream(end2)
    
    gotOutput = []
    for r1,r2 in BEDDuplicateFilter(end1_infs, end2_infs, removeJuncTags=True, removePETags=True):
      r1 = [x.name for x in r1]
      r2 = [x.name for x in r2]
      gotOutput.append((r1, r2))
    
    if debug :
      print "got \n------ "
      for r1,r2 in gotOutput : print str(r1) + "," + str(r2)
      print "-------"
      print "expect \n------ "
      for r1,r2 in expectOutput : print str(r1) + "," + str(r2)
      print "-------"
    
    expectOutput.sort()                                    
    gotOutput.sort()
    self.assertTrue(gotOutput == expectOutput)
    
  def testBEDSortedByReadName(self):
    sortedIn = "chr1 \t 32 \t 43 \t readA\n" +\
               "chr2 \t 23 \t 57 \t readB\n" +\
               "chr1 \t 32 \t 46 \t readC\n" +\
               "chr5 \t 76 \t 99 \t readD\n" +\
               "chr5 \t 32 \t 45 \t readE\n" +\
               "chr1 \t 91 \t 97 \t readF\n" +\
               "chr6 \t 32 \t 43 \t readG\n" +\
               "chr7 \t 53 \t 83 \t readH\n"
    unsortedIn = "chr1 \t 32 \t 43 \t readF\n" +\
                 "chr2 \t 23 \t 57 \t readB\n" +\
                 "chr3 \t 32 \t 46 \t readA\n" +\
                 "chr1 \t 76 \t 99 \t readD\n" +\
                 "chr5 \t 32 \t 45 \t readE\n" +\
                 "chr5 \t 91 \t 97 \t readM\n" +\
                 "chr1 \t 32 \t 43 \t readC\n" +\
                 "chr6 \t 53 \t 83 \t readH\n"
    sorteds = DummyInputStream(sortedIn)
    unsorteds = DummyInputStream(unsortedIn)
    
    def run(strm):
      for item in BEDIterator(strm, sortedby = ITERATOR_SORTED_NAME) : pass
      
    # should run without any problems...
    run(sorteds)
    # should raise an exception..
    self.assertRaises(BEDError, run, unsorteds)
    
  def testBEDSortedByChrom(self):
    sortedIn = "chr1 \t 32 \t 43 \t readB\n" +\
               "chr1 \t 23 \t 57 \t readA\n" +\
               "chr1 \t 32 \t 46 \t readF\n" +\
               "chr5 \t 76 \t 99 \t readH\n" +\
               "chr5 \t 32 \t 45 \t readE\n" +\
               "chr6 \t 91 \t 97 \t readC\n" +\
               "chr6 \t 32 \t 43 \t readG\n" +\
               "chr6 \t 53 \t 83 \t readD\n"
    unsortedIn = "chr1 \t 32 \t 43 \t readA\n" +\
                 "chr2 \t 23 \t 57 \t readB\n" +\
                 "chr3 \t 32 \t 46 \t readC\n" +\
                 "chr1 \t 76 \t 99 \t readD\n" +\
                 "chr5 \t 32 \t 45 \t readE\n" +\
                 "chr5 \t 91 \t 97 \t readF\n" +\
                 "chr1 \t 32 \t 43 \t readG\n" +\
                 "chr6 \t 53 \t 83 \t readH\n"
    sorteds = DummyInputStream(sortedIn)
    unsorteds = DummyInputStream(unsortedIn)
    
    def run(strm):
      for item in BEDIterator(strm, sortedby = ITERATOR_SORTED_CHROM) : pass
      
    # should run without any problems...
    run(sorteds)
    # should raise an exception..
    self.assertRaises(BEDError, run, unsorteds)
    
  def testBEDSortedByStart(self):
    sortedIn = "chr1 \t 32 \t 43 \t readF\n" +\
               "chr1 \t 43 \t 57 \t readB\n" +\
               "chr1 \t 44 \t 46 \t readH\n" +\
               "chr5 \t 23 \t 30 \t readD\n" +\
               "chr5 \t 30 \t 45 \t readE\n" +\
               "chr6 \t 22 \t 30 \t readA\n" +\
               "chr6 \t 23 \t 29 \t readG\n" +\
               "chr7 \t 53 \t 83 \t readC\n"
    unsortedIn = "chr1 \t 32 \t 43 \t readF\n" +\
                 "chr1 \t 30 \t 57 \t readB\n" +\
                 "chr1 \t 44 \t 46 \t readH\n" +\
                 "chr5 \t 23 \t 30 \t readD\n" +\
                 "chr8 \t 30 \t 45 \t readE\n" +\
                 "chr6 \t 22 \t 30 \t readA\n" +\
                 "chr6 \t 10 \t 29 \t readG\n" +\
                 "chr7 \t 53 \t 83 \t readC\n"
    sorteds = DummyInputStream(sortedIn)
    unsorteds = DummyInputStream(unsortedIn)
    
    def run(strm):
      for item in BEDIterator(strm, sortedby = ITERATOR_SORTED_START) : pass
      
    # should run without any problems...
    run(sorteds)
    # should raise an exception..
    self.assertRaises(BEDError, run, unsorteds)
    
      
  def testSizeOfOverlap(self):
    MAX_INDEX = 100000
    ITERATIONS = 100
    
    for i in range(0, ITERATIONS) :
      region1 = randomBEDElement(maxIndex = MAX_INDEX)
      region2 = copy.deepcopy(region1)
      
      ## region 2 entirely inside region 1 -- ans = size of region 2
      region2.start = int(random.random() * (region1.end - region1.start - 1) +\
                      region1.start)
      region2.end = int(random.random() * (region1.end - region2.start) +\
                    region2.start)
      self.assertEqual(region1.sizeOfOverlap(region2), len(region2))
      self.assertEqual(region1.sizeOfOverlap(region2), 
                       region2.sizeOfOverlap(region1))
      
      ## region 2 start inside region 1, but end beyond it -- ans = region1.end - region2.start
      region2.start = int(random.random() * (region1.end - region1.start - 1) +\
                      region1.start)
      region2.end = int(random.random() * (MAX_INDEX - region1.end) +\
                    region1.end)
      self.assertEqual(region1.sizeOfOverlap(region2), region1.end - region2.start + 1)
      self.assertEqual(region1.sizeOfOverlap(region2), 
                       region2.sizeOfOverlap(region1))
      
      ## region 2 end inside region 1, but start beyond it -- ans = region2.end - region1.start
      region2.end = int(random.random() * (region1.end - region1.start - 1) +\
                    region1.start)
      region2.start = int(random.random() * (region1.start))
      self.assertEqual(region1.sizeOfOverlap(region2), region2.end - region1.start + 1)
      self.assertEqual(region1.sizeOfOverlap(region2), 
                       region2.sizeOfOverlap(region1))
      
      ## no intersection at all -- ans = 0
      region2.start = int(random.random() * (region1.start - 1))
      region2.end = int(random.random() * (region1.start - region2.start) +\
                    region2.start)
      self.assertEqual(region1.sizeOfOverlap(region2), 0)
      self.assertEqual(region1.sizeOfOverlap(region2), 
                       region2.sizeOfOverlap(region1))
  
  def testSubtract(self):
    pass
  
  def testDistance(self):
    pass
  

if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])
