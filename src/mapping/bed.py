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
  History:       
  
"""  
"""              18th August 2010 -- Philip Uren
                   * added BEDIterator 
                   * updated header details 
                   * added testSizeOfOverlap unit test
"""
"""              20th August 2010 -- Philip Uren
                   * added check for same chrom to intersection test
"""
"""              13th September 2010 -- Philip Uren
                   * modified BEDElement and BEDElementFromString to allow 
                     BED elements to have only chrom, start and end 
                   * cleaned details of junction reads
"""
"""              16th September 2010 -- Philip Uren
                   * added check for sorted order to BEDIterator
""" 
"""              17th September 2010 -- Philip Uren
                   * added intervalTree function
""" 
"""              24th September 2010 -- Philip Uren
                   * fixed bugs with checking for sorted BED files when 
                     not required
"""
"""              1st October 2010 -- Philip Uren
                   * added ability to check BED file is sorted by read name
                   * added BEDUniqueFilter
                   * added unit tests for above two items
                   * added missing exception class 
                   * added unit test for BED file being sorted by start
                   * added verbose option to BEDIterator
"""
"""              2nd October 2010 -- Philip Uren
                   * added BEDDuplicateFilter and associated unit tests
"""
"""              19th October 2010 -- Philip Uren
                   * added best option for BEDUniqueFilter
                   * changed type for score from float to int
"""                    
"""              22nd October 2010 -- Philip Uren
                   * fixed bug in BEDIterator when passing string instead of
                     stream 
                   * added detail to the exception raised when BED elements 
                     appear not to have enough parts to them
""" 
"""              23rd October 2010 -- Philip Uren
                   * added BEDIterator sorting order of chromosome 
                   * added unit test for above
""" 
"""              27th October 2010 -- Philip Uren
                   * moved to smithlab_py 
                   * changed usage of progress indicator to use new class
"""
"""              28th October 2010 -- Philip Uren 
                   * changed unit test for overlap
                   * removed iterators to a new module
""" 
"""              16th November 2010 -- Philip Uren
                   * added support for colours and more error checking on 
                     parsing. Cleaned up some old code
"""
"""              24th Novemeber 2010 -- Philip Uren
                   * added toGenomicCoordinates and associated unit tests
"""
"""              20th December 2010 -- Philip Uren
                   * allow parsing to ignore data that doesn't match BED spec 
                     (e.g. things that should be numbers, but aren't)
"""
"""              06th January 2011 -- Philip Uren
                   * moved intervalTree method from this file to 
                     mapping.bedIterators
"""
"""              28th December 2011 -- Philip Uren
                   * simplified element subtraction code, allowed multiple 
                     elements to be subtracted at same time, added unit tests
                     for subtraction
""" 

"""  
  TODO:         
                * finish unit tests
                * lots of duplication between intervalTreesFromList and 
                  intervalTrees that should be removed  
"""

import sys, os, unittest, copy

from testing.dummyfiles import DummyInputStream, DummyOutputStream
from util.progressIndicator import ProgressIndicator
from util.fileUtils import linesInFile
from datastruct.intervalTree import IntervalTree

class BEDError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)

def intervalTreesFromList(inElements, verbose = False):
  elements = {}
  if verbose :
    totalLines = len(inElements)
    pind = ProgressIndicator(totalToDo = totalLines, 
                                   messagePrefix = "completed", 
                                   messageSuffix = "of parsing")      
  
  for element in inElements :
    if not element.chrom in elements : elements[element.chrom] = []
    elements[element.chrom].append(element)
    if verbose :
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


def toGenomicCoordinates(start, end, transcript, debug = False):
  """
    @summary: transform transcript coordinates into genomic coordinates,
              returning a BEDElement (or two if the region spans a splice
              junction)
    @param start: the start location in transcript coordinates
    @param end: the end location in transcript coordinates
    @param transcript: list of BEDElement objects describing where
                       the exons of this transcript are
    @return: a list of BEDElement objects representing the regions of the 
             transcript covered in genomic coordinates -- this may be only 
             one region if it's inside a single exon, or multiple if spanning
             more than one exon
    @raise BEDError: if the transcript has exons on multiple chromosomes
  """
  
  def regionIterator (start, end, transcript, debug = False):
    """
      @summary: an internal function for iterating covered regions
    """
    # do some quick error checking..
    totalSize = sum([x.end - x.start for x in transcript])
    if start < 0 or end > totalSize :
      raise BEDError("gave start and end as " + str(start) + "," +str(end) +\
                     " -- transcript is " + str(totalSize) + " long")
    
    # sort the exons into order and figure out where the gene is
    transcript.sort(key = lambda x : x.start)
    chrom = None 
    for exon in transcript : 
      if chrom == None : chrom = exon.chrom
      elif chrom != exon.chrom : 
        raise BEDError ("got transcript with exons on multiple chromosomes!")
    tStart = transcript[0].start
    
    # figure out the relative location of the splice junctions
    cumulative = 0
    junctions = [0]
    for exon in transcript :
      cumulative += exon.end - exon.start
      junctions.append(cumulative)
      
    if debug :
      sys.stderr.write("identified exon boundaries at: "+str(junctions)+"\n")
      
    # now figure out which exon the start and end are in
    startInd, endInd = None, None
    prev = 0
    exonInd = 0
    for i in range(1,len(junctions)) :
      junc = junctions[i]
      
      # the start...
      if debug : sys.stderr.write("is " + str(start) + " in " + str(prev) +\
                                  "-" + str(junc) + "? ")
      if start >= prev and start < junc :
        startInd = exonInd
        if debug : sys.stderr.write("yes \n")
      elif debug : sys.stderr.write("no \n")
      
      # the end...
      if debug : sys.stderr.write("is " + str(end) + " in " + str(prev) +\
                                  "-" + str(junc) + "? ")
      if end >= prev and end <= junc :
        endInd = exonInd
        if debug : sys.stderr.write("yes \n")
      elif debug : sys.stderr.write("no \n")
      
      exonInd += 1
      prev = junc
    
    # convert to absolute coordinates
    absStart = (start - junctions[startInd]) + transcript[startInd].start
    absEnd = (end - junctions[endInd]) + transcript[endInd].start
    
    # first we yeild any exon that is entirely covered by the region
    for exon in transcript :
      if absStart < exon.start and absEnd > exon.end :
        yield BEDElement(exon.chrom, exon.start, exon.end)
    
    # if both start and end are in same exon, we can yield a single region
    if startInd == endInd :
      yield BEDElement(chrom, absStart, absEnd) 
    else :
      yield BEDElement(chrom, absStart, transcript[startInd].end)
      yield BEDElement(chrom, transcript[endInd].start, absEnd)
  
  return [region for region in regionIterator(start, end, transcript, debug)]
          

def BEDElementFromString(line, scoreType=int, dropAfter = None):
  """
    @summary: Given a string in BED format, parse the string and return a 
              BEDElement object
    @param line: the string to be parsed
    @param dropAfter: an int indicating that any fields after and including this
                      field should be ignored as they don't conform to the BED 
                      format. By default, None, meaning we use all fields. Index
                      from zero.
    @return: BEDElement 
  """
  peices = line.split("\t")
  if dropAfter != None : peices = peices[0:dropAfter]
  if len(peices) < 3 : 
    raise BEDError("BED elements must have at least chrom, start and end" +\
                   " found only " + str(len(peices)) + " in " + line)
  chrom = peices[0]
  start = peices[1]
  end = peices[2]
  
  name = None 
  score = None 
  strand = None 
  thickStart = None
  thickEnd = None 
  colour = None 
  blockCount = None 
  blockSizes = None 
  blockStarts = None 
  
  if len(peices) >= 4 != None : name = peices[3]
  if len(peices) >= 5 != None : score = peices[4]
  if len(peices) >= 6 != None : strand = peices[5]
  if len(peices) >= 7 != None : thickStart = peices[6]
  if len(peices) >= 8 != None : thickEnd = peices[7]
  if len(peices) >= 9 != None : colour = peices[8]
  if len(peices) >= 10 != None : blockCount = peices[9]
  if len(peices) >= 11 != None : blockSizes = peices[10]
  if len(peices) >= 12 != None : blockStarts = peices[11]
  
  # some programs output things that are meant to be BED format, 
  # but don't obey the rules... we're going to ignore thickStart 
  # and thickEnd values that are not ints
  try : int(thickStart)
  except : thickStart = None 
  try : int(thickEnd)
  except : thickEnd = None 
  
  return BEDElement(chrom, start, end, name, score, strand, thickStart,
                    thickEnd, colour, blockCount, blockSizes, blockStarts, scoreType)
  

class BEDElement :
  POSITIVE_STRAND = "+"
  NEGATIVE_STRAND = "-"
  DEFAULT_STRAND = POSITIVE_STRAND
  
  def __init__(self, chrom, start, end, name = None, score = None, 
               strand = None, thickStart = None, thickEnd = None, 
               colour = None, blockCount = None, blockSizes = None,
               blockStarts = None, scoreType = int):
    """
      @summary: Constructor for the BEDElement class
      @note: only the first three parameters are required 
      @note: if any parameter is omitted, no default will be set for it
             (it will just be = None) and it won't appear in any output.
      @note: if any parameter is provided, all parameters that proceed 
             it must also be provided. 
      @note: BEDElements are inclusive of the start, but not the end 
             coordinate 
    """
    # the basic read info -- we must get at least this much 
    if chrom == None or start == None or end == None :
      raise BEDError("Must provided at least chrom, start, end for BED element")
    self.chrom = chrom.strip()
    self.start = int(start)
    self.end = int(end)
    
    # we might get the following too...
    # name:
    self.name = None
    if name != None: self.name = name.strip()
    
    # score
    if score != None and name == None :
      raise BEDError("If score is provided, must also provide name")
    self.score = None
    if score != None : self.score = scoreType(score)
    
    # strand 
    if strand != None and (name == None or score == None) :
      raise BEDError("If strand is provided, must also provide name and score")
    self.strand = None  
    if strand != None : self.strand = strand.strip()
    
    # thickStart
    if thickStart != None and (name == None or score == None or strand == None) :
      raise BEDError("If thickStart is provided, must also provide name, " +\
                     "score and strand")
    self.thickStart = None
    if thickStart != None : self.thickStart = int(thickStart)
    
    # thickEnd
    if thickEnd != None and (name == None or score == None or strand == None \
                             or thickStart == None) :
      raise BEDError("If thickEnd is provided, must also provide name, " +\
                     "score, strand and thickStart")
    self.thickEnd = None
    if thickEnd != None : self.thickEnd = int(thickEnd)
    
    # colour
    if colour != None and (name == None or score == None or strand == None \
                             or thickStart == None or thickEnd == None) :
      raise BEDError("If colour is provided, must also provide name, " +\
                     "score, strand, thickStart and thickEnd")
    self.colour = None
    if colour != None : self.colour = colour
    
    # blockCount
    if blockCount != None and (name == None or score == None or strand == None \
                             or thickStart == None or thickEnd == None \
                             or colour == None) :
      raise BEDError("If colour is provided, must also provide name, " +\
                     "score, strand, thickStart, thickEnd, colour")
    self.blockCount = None
    if blockCount != None : self.blockCount = blockCount
    
    # blockSizes
    if blockSizes != None and (name == None or score == None or strand == None \
                             or thickStart == None or thickEnd == None \
                             or colour == None or blockCount == None) :
      raise BEDError("If blockSizes is provided, must also provide name, " +\
                     "score, strand, thickStart, thickEnd, colour and " +\
                     "blockCount")
    self.blockSizes = None
    if blockSizes != None : self.blockSizes = blockSizes
    
    # blockStarts
    if blockStarts != None and (name == None or score == None or strand == None \
                             or thickStart == None or thickEnd == None \
                             or colour == None or blockCount == None \
                             or blockSizes == None) :
      raise BEDError("If blockStarts is provided, must also provide name, " +\
                     "score, strand, thickStart, thickEnd, colour, " +\
                     "blockCount and blockSizes")
    self.blockStarts = None
    if blockStarts != None : self.blockStarts = blockStarts

  
  def __eq__(self, e):
    if e == None : return False
    return self.chrom == e.chrom and self.start == e.start and\
           self.end == e.end and self.name == e.name and\
           self.score == e.score and self.strand == e.strand
           
  def __lt__(self, rhs):
    """
      @summary: is self < rhs? default comparison by end 
    """
    if rhs == None : return False
    if self.chrom < rhs.chrom : return True
    if self.chrom > rhs.chrom : return False
    if self.end < rhs.end : return True
    return False
           
  def sameRegion(self, e):
    """
      @summary: return true if self and e are for the same region 
                (ignores differences in non-region related fields)
    """
    if e == None : return False
    return self.chrom == e.chrom and self.start == e.start and\
           self.end == e.end and self.name == e.name and\
           self.strand == e.strand
    
        
  def __len__(self):
    return (self.end - self.start)       
        
  def __str__(self):
    """
      @summary: Produce a string representation of the BED element. Only 
                those fields which we have values for will be output. 
    """
    delim = "\t"
    res = self.chrom + delim + str(self.start) + delim + str(self.end)
    if self.name != None : res += (delim + str(self.name))
    if self.score != None : res += (delim + str(self.score))
    if self.strand != None : res += (delim + str(self.strand))
    if self.thickStart != None : res += (delim + str(self.thickStart))
    if self.thickEnd != None : res += (delim + str(self.thickEnd))
    if self.colour != None : res += (delim + str(self.colour))
    if self.blockCount != None : res += (delim + str(self.blockCount))
    if self.blockSizes != None : res += (delim + str(self.blockSizes))
    if self.blockStarts != None : res += (delim + str(self.blockStarts))
    
    return  res
  
  def distance(self, e):
    dist = 0
    if not e.intersects(self) :
      dist = max(self.start, e.start) - min(self.end, e.end)
    return dist
  
  
  def __singleIntersect(self, e):
    """
      @summary: this is a private method which will return a list of regions
                that represent the result of subtracting e from self. The 
                list might be empty if self is totally contained in e, or 
                may contain two elements if e is totally contained in self. 
                otherwise there'll be one element.
    """
    if e.chrom != self.chrom or e.end < self.start or e.start > self.end :
      # no intersection
      return [copy.copy(self)]
    if e.start <= self.start and e.end >= self.end :
      # whole self is removed.
      return []
    if e.start > self.start and e.end < self.end :
      # splits self in two
      r1 = copy.copy(self)
      r2 = copy.copy(self)
      r1.end = e.start
      r2.start = e.end
      return [r1,r2]
    if e.start <= self.start : 
      # cuts off start of self
      r = copy.copy(self)
      r.start = e.end
      return [r]
    if e.end >= self.end :
      # cuts off end of self
      r = copy.copy(self)
      r.end = e.start
      return [r]
    # oops, we screwed up if we got to here.. 
    raise BEDError("fatal error - failed BED subtraction of " +\
                   str(e) + " from " + str(self))
    
  
  def subtract(self, es):
    """ 
      @summary: subtracts the BED elements in es from self
      @param es: a list of BED elements (or anything with chrom, start, end)
      @return: a list of BED elements which represent what is left of 
               self after the subtraction. This might be an empty list.
    """
    workingSet = [self]
    for e in es : 
      newWorkingSet = []
      for w in workingSet :
        newWorkingSet += w.__singleIntersect(e)
      workingSet = newWorkingSet
    return workingSet

    
  def sizeOfOverlap(self, e):
    # no overlap
    if not self.intersects(e) : return 0
    
    # complete inclusion..
    if e.start >= self.start and e.end <= self.end : return len(e)
    if self.start >= e.start and self.end <= e.end : return len(self)
    
    # partial overlap
    if e.start > self.start : return (self.end - e.start)
    if self.start > e.start : return (e.end - self.start)
    
    
  def intersects(self, e):
    """ returns true if this elements intersects the element e """

    if self.chrom != e.chrom : return False
    if self.end >= e.start and self.end <= e.end : return True
    if self.start >= e.start and self.start <= e.end : return True
    if e.start >= self.start and e.start <= self.end : return True
    if e.end >= self.start and e.end <= self.end : return True
    return False
  
  def isPositiveStrand(self):
    if self.strand == None and self.DEFAULT_STRAND == self.POSITIVE_STRAND :
      return True
    return self.strand == self.POSITIVE_STRAND 
  
  def isNegativeStrand(self):
    if self.strand == None and self.DEFAULT_STRAND == self.NEGATIVE_STRAND :
      return True
    return self.strand == self.NEGATIVE_STRAND
  
  def transcriptKey(self):
    """
      @summary: returns a string that uniquely identifies the name, chrom and
                strand of this element, assuming this identifies an exon and 
                has a name of the format refSeq_exon_details
    """
    refseq = self.name
    if self.name.find("_exon_") != -1 : refseq = self.name.split("_exon_")
    return refseq + self.strand + self.chrom
     

class BEDUnitTests(unittest.TestCase):
  """
    @summary: Unit tests for functions and classes in this module  
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
    debug = False
    a = BEDElement("chr1",10,100)
    b1 = BEDElement("chr1",1,8)
    b2 = BEDElement("chr1",15,20)
    b3 = BEDElement("chr1",50,60)
    b4 = BEDElement("chr1",90,120)
    c1 = BEDElement("chr2",30,40)
    
    res = a.subtract([b1,b2,b3,b4,c1])
    expect = [BEDElement("chr1",10,15),
              BEDElement("chr1",20,50),
              BEDElement("chr1",60,90)]
    res.sort()
    expect.sort()
    if debug :
      sys.stderr.write("\nDebugging BED subtraction test\n")
      sys.stderr.write("expect\n" + "\n".join([str(x) for x in expect]) + "\n")
      sys.stderr.write("got\n" + "\n".join([str(x) for x in res]) + "\n")
    assert(res == expect)
  
  def testDistance(self):
    pass
  
  def testToGenomicCoordinates(self):
    debug = False
    exon1 = BEDElement("chr1", 101, 110)
    exon2 = BEDElement("chr1", 120, 130)
    exon3 = BEDElement("chr1", 170, 175)
    
    transcript = [exon1,exon2,exon3]
    sortKey = lambda x : x.start
    
    res1 = toGenomicCoordinates(8, 15, transcript, debug = debug)
    res2 = toGenomicCoordinates(8, 24, transcript, debug = debug)
    res3 = toGenomicCoordinates(0, 3, transcript, debug = debug)
    for res in [res1,res2,res3] : res.sort(key = sortKey)
    
    expect1 = [BEDElement("chr1", 109, 110), BEDElement("chr1", 120, 126)]
    expect2 = [BEDElement("chr1", 109, 110), BEDElement("chr1", 120, 130), 
               BEDElement("chr1", 170, 175)]
    expect3 = [BEDElement("chr1", 101, 104)]
    for expect in [expect1, expect2, expect3] : expect.sort(key = sortKey)
    
    if debug :
      for expect, got in [(expect1, res1), (expect2, res2), (expect3, res3)] :
        sys.stderr.write("expect: \n")
        for e in expect : sys.stderr.write("\t" + str(e) + "\n")
        sys.stderr.write("got: \n")
        for e in got : sys.stderr.write("\t" + str(e) + "\n")
        sys.stderr.write("\n")
    
    self.assertTrue(res1 == expect1)
    self.assertTrue(res2 == expect2)
    self.assertTrue(res3 == expect3)

if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])
