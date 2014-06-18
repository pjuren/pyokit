#!/usr/bin/python

"""
  Date of Creation: 3rd Apr 2010    
                       
  Description:   Classes and functions for manipulating Genomic Intervals 

  Copyright (C) 2010  
  University of Southern California,
  Philip J. Uren,
  
  Authors: Philip J. Uren
  
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
"""
  

import sys, os, unittest, copy

from pyokit.testing.dummyfiles import DummyInputStream, DummyOutputStream
from pyokit.util.progressIndicator import ProgressIndicator
from pyokit.util.fileUtils import linesInFile
from pyokit.datastruct.intervalTree import IntervalTree


###############################################################################
##                           EXCEPTION CLASSES                               ##
###############################################################################

class GenomicIntervalError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)
  
  
###############################################################################
##       FUNCTIONS FOR MANIPULATING COLLECTIONS OF GENOMIC INTERVALS         ##
###############################################################################

def intervalTreesFromList(inElements, verbose = False):
  """
    @summary: build a dictionary, indexed by chromosome name, of interval trees
              for each chromosome.
    @param inElements: list of genomic intervals. Members of the list must have
                       chrom, start and end fields; no other restrictions. 
    @param verbose: output progress messages to sys.stderr if True 
  """
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


def collapseRegions(s):
  """
    @summary: given a list of genomic intervals with chromosome, start and end  
              field, collapse into a set of non-overlapping intervals. Intervals
              must be sorted by chromosome and then start coordinate.
    @return:  list of intervals that define the collapsed regions. Note that
              these are all new objects, no existing object from s is returned
              or altered. Returned regions will all have name "X", strand +
              and score 0
    @param s: list of genomic regions to collapse 
    @raise BEDError: if the input regions are not correctly sorted (chromosome 
                     then start) 
    @note: O(n) time, O(n) space
  """
  debug = False
  
  if len(s) == 0 or len(s) == 1 : return copy.deepcopy(s)
  
  res = []
  current = copy.copy(s[0])
  current.strand = '+'
  current.score = 0
  current.name = "X"
  for i in range(1, len(s)) :
    if debug :
      sys.stderr.write("processing " + str(s[i]) + "\n")
      sys.stderr.write("\tcurrent is: " + str(current) + "\n")
    
    # make sure things are sorted.. 
    if (s[i].chrom < s[i-1].chrom) or \
       (s[i].chrom == s[i-1].chrom and s[i].start < s[i-1].start) :
      raise BEDError("collapsing regions failed. saw this region: " +\
                     str(s[i-1]) + " before this one: " + str(s[i]))
      
    # because of sorting order, we know that nothing else exists with 
    # start less than s[i] which we haven't already seen.
    if s[i].start > current.end or s[i].chrom != current.chrom : 
      res.append(current)
      current = copy.copy(s[i])
      current.strand = '+'
      current.score = 0
      current.name = "X"
    else :
      current.end = max(s[i].end, current.end)
      
  # don't forget the last one...
  res.append(current)
  
  return res
    
    
def regionsIntersection(s1, s2):
  """
    @summary: given two lists of genomic regions with chromosome, start and end 
              coordinates, return a new list of regions which is the 
              intersection of those two sets. Lists must be sorted by 
              chromosome and start index
    @return: new list that represents the intersection of the two input lists.
             output regions will all have name "X", be one strand "+" and have
             score 0
    @param s1: first list of genomic regions
    @param s2: second list of genomic regions
    @raise BEDError: if the input regions are not sorted correctly (by 
                     chromosome and start index)
    @note: O(n) time, O(n) space; informally, might use up to 3x space of 
           input 
  """
  debug = False
  
  # we don't need to explicitly check for sorting because sorted order is 
  # a post-condition of the collapsing function
  s1_c = collapseRegions(s1)
  s2_c = collapseRegions(s2)
  
  if len(s1_c) == 0 or len(s2_c) == 0 : return []
  
  res = []
  j = 0
  for i in range(0, len(s1_c)) :
    if debug :
      sys.stderr.write("processing from s1_c : " + str(s1_c[i]) + "\n")
    
    # find first thing in s2_c with end in or after s1_c[i]
    hits = True
    if debug : sys.stderr.write("i = " + str(i) + " and j = " + str(j) + "\n")
    while j < len(s2_c) and \
          (s2_c[j].chrom < s1_c[i].chrom or \
          (s2_c[j].chrom == s1_c[i].chrom and s2_c[j].end <= s1_c[i].start)) : 
      j += 1
    # nothing intersects if we hit the end of s2, or the end of the chrom, 
    # or we're still on the same chrom but start after the end of s2_c[i]
    if j >= len(s2_c) or s2_c[j].chrom > s1_c[i].chrom or \
       (s2_c[j].chrom == s1_c[i].chrom and s2_c[j].start >= s1_c[i].end) : 
      continue 
    
    # now everything at or after j in s2_c that starts before
    # the end of s1_c must overlap with it
    while s2_c[j].start < s1_c[i].end :
      s = max(s1_c[i].start, s2_c[j].start)
      e = min(s1_c[i].end, s2_c[j].end)
      overlap = GenomicInterval(s1_c[i].chrom,s,e,"X",0,"+")
      if debug : sys.stderr.write("\tadding to overlaps: " +\
                                  str(overlap) + "\n")
      res.append(overlap)
      j += 1
      if j >= len(s2_c) or s2_c[j].chrom != s1_c[i].chrom : break
    
    # it's possible the last intersecting element runs on to the 
    # next element from s1_c, so...
    j -= 1
    if debug : sys.stderr.write("\tmoving s2_c index back to " +\
                                str(s2_c[j]) + "\n")
  
  return res
          
          
###############################################################################
##                FUNCTIONS FOR PARSING GENOMIC INTERVALS                    ##
###############################################################################

def parseBEDString(line, scoreType=int, dropAfter = None):
  """
    @summary: Given a string in BED format, parse the string and return a 
              GenomicInterval object
    @param line: the string to be parsed
    @param dropAfter: an int indicating that any fields after and including this
                      field should be ignored as they don't conform to the BED 
                      format. By default, None, meaning we use all fields. Index
                      from zero.
    @return: GenomicInterval object built from the BED string representation 
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
  
  return GenomicInterval(chrom, start, end, name, score, strand, thickStart,
                         thickEnd, colour, blockCount, blockSizes, blockStarts, 
                         scoreType)


###############################################################################
##                          GENOMIC INTERVAL CLASS                           ##
###############################################################################

class GenomicInterval :
  POSITIVE_STRAND = "+"
  NEGATIVE_STRAND = "-"
  DEFAULT_STRAND = POSITIVE_STRAND
  
  def __init__(self, chrom, start, end, name = None, score = None, 
               strand = None, thickStart = None, thickEnd = None, 
               colour = None, blockCount = None, blockSizes = None,
               blockStarts = None, scoreType = int):
    """
      @summary: Constructor for the GenomicInterval class
      @note: only the first three parameters are required 
      @note: if any parameter is omitted, no default will be set for it
             (it will just be = None) and it won't appear in any output.
      @note: if any parameter is provided, all parameters that proceed 
             it must also be provided. 
      @note: GenomicIntervals are inclusive of the start, but not the end 
             coordinate 
    """
    # the basic read info -- we must get at least this much 
    if chrom == None or start == None or end == None :
      raise GenomicIntervalError("Must provided at least chrom, start, end " +\
                                 "for Genomic Interval")
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
      raise GenomicIntervalError("If blockStarts is provided, must also "    +\
                                 "provide name, score, strand, thickStart, " +\
                                 "thickEnd, colour, blockCount and blockSizes")
    self.blockStarts = None
    if blockStarts != None : self.blockStarts = blockStarts

  
  def __hash__(self):
    """
      @summary: return a hash of this GenomicInterval
    """
    return hash(str(self))
  
  def __eq__(self, e):
    """
      @summary: return true if two GenomicInterval objects are equal
    """
    if e == None : return False
    try :
      return  self.chrom == e.chrom and self.start == e.start and\
              self.end == e.end and self.name == e.name and\
              self.score == e.score and self.strand == e.strand
    except AttributeError :
      return False
           
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
    """
      @summary: return the distance from this GenomicInterval to e. We consider
                intervals that overlap to have a distance of 0 to each other. 
                The distance between two intervals on different chromosomes is 
                considered undefined, and causes an exception to be raised.
      
    """
    if e.chrom != self.chrom :
      raise GenomicIntervalError("cannot get distance from " + str(self) +\
                                 " to " + str(e) + " as they're on "     +\
                                 "different chromosomes")
    dist = 0
    if not e.intersects(self) :
      dist = max(self.start, e.start) - min(self.end, e.end)
    return dist
  
  def signedDistance(self, e):
    """ 
      @summary: return the distance from self to e, with sign. If e comes 
                earlier than self, the distance will be negative. We consider
                intervals that overlap to have a distance of 0 to each other. 
                The distance between two intervals on different chromosomes is 
                considered undefined, and causes an exception to be raised.
    """
    dist = self.distance(e)
    if e < self : dist = dist * -1
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
    """
      @summary: returns the number of bases that are shared in common 
                between self and e.
    """
    # no overlap
    if not self.intersects(e) : return 0
    
    # complete inclusion..
    if e.start >= self.start and e.end <= self.end : return len(e)
    if self.start >= e.start and self.end <= e.end : return len(self)
    
    # partial overlap
    if e.start > self.start : return (self.end - e.start)
    if self.start > e.start : return (e.end - self.start)
    
    
  def intersects(self, e):
    """ 
      @summary: returns true if this elements intersects the element e 
    """

    if self.chrom != e.chrom : return False
    if self.end >= e.start and self.end <= e.end : return True
    if self.start >= e.start and self.start <= e.end : return True
    if e.start >= self.start and e.start <= self.end : return True
    if e.end >= self.start and e.end <= self.end : return True
    return False
  
  def isPositiveStrand(self):
    """
      @summary: returns true if this element is on the positive strand
    """
    if self.strand == None and self.DEFAULT_STRAND == self.POSITIVE_STRAND :
      return True
    return self.strand == self.POSITIVE_STRAND 
  
  def isNegativeStrand(self):
    """
      @summary: returns true if this element is on the negative strand
    """
    if self.strand == None and self.DEFAULT_STRAND == self.NEGATIVE_STRAND :
      return True
    return self.strand == self.NEGATIVE_STRAND
  
  
###############################################################################
##                        UNIT TESTS FOR THIS MODULE                         ##
###############################################################################

class GenomicIntervalUnitTests(unittest.TestCase):
  """
    @summary: Unit tests for functions and classes in this module  
  """
  
  def setUp(self):
    pass
  
  def testCollapse(self):
    debug = False
    elements = ["chr1"+"\t"+"10"+"\t"+"20"+"\t"+"R01"+"\t"+"0"+"\t"+"+",
                "chr1"+"\t"+"30"+"\t"+"40"+"\t"+"R02"+"\t"+"1"+"\t"+"+",
                "chr1"+"\t"+"35"+"\t"+"50"+"\t"+"R03"+"\t"+"0"+"\t"+"+",
                "chr1"+"\t"+"45"+"\t"+"65"+"\t"+"R04"+"\t"+"0"+"\t"+"+",
                "chr1"+"\t"+"55"+"\t"+"60"+"\t"+"R05"+"\t"+"3"+"\t"+"-",
                "chr1"+"\t"+"70"+"\t"+"80"+"\t"+"R06"+"\t"+"0"+"\t"+"+",
                "chr1"+"\t"+"75"+"\t"+"95"+"\t"+"R07"+"\t"+"0"+"\t"+"+",
                "chr1"+"\t"+"85"+"\t"+"90"+"\t"+"R08"+"\t"+"1"+"\t"+"-",
                "chr2"+"\t"+"40"+"\t"+"60"+"\t"+"R10"+"\t"+"0"+"\t"+"+",
                "chr3"+"\t"+"10"+"\t"+"20"+"\t"+"R11"+"\t"+"4"+"\t"+"+",
                "chr3"+"\t"+"20"+"\t"+"30"+"\t"+"R12"+"\t"+"0"+"\t"+"-"]
    expect_str = ["chr1"+"\t"+"10"+"\t"+"20"+"\t"+"X"+"\t"+"0"+"\t"+"+",
                  "chr1"+"\t"+"30"+"\t"+"65"+"\t"+"X"+"\t"+"0"+"\t"+"+",
                  "chr1"+"\t"+"70"+"\t"+"95"+"\t"+"X"+"\t"+"0"+"\t"+"+",
                  "chr2"+"\t"+"40"+"\t"+"60"+"\t"+"X"+"\t"+"0"+"\t"+"+",
                  "chr3"+"\t"+"10"+"\t"+"30"+"\t"+"X"+"\t"+"0"+"\t"+"+"]
    input = [parseBEDString(x) for x in elements]
    expect = [parseBEDString(x) for x in expect_str]
    got = collapseRegions(input)
    if debug :
      sys.stderr.write("expect:\n")
      sys.stderr.write("\n".join([str(x) for x in expect]) + "\n")
      sys.stderr.write("got:\n")
      sys.stderr.write("\n".join([str(x) for x in got]) + "\n")
    self.assertEqual(expect, got)
    
  def testRegionsIntersection(self):
    debug = False
    s1_elements = ["chr1"+"\t"+"40" +"\t"+"90" +"\t"+"R11" +"\t"+"0"+"\t"+"+",
                   "chr1"+"\t"+"100"+"\t"+"120"+"\t"+"R12" +"\t"+"3"+"\t"+"+",
                   "chr1"+"\t"+"160"+"\t"+"190"+"\t"+"R13" +"\t"+"1"+"\t"+"-",
                   "chr1"+"\t"+"200"+"\t"+"210"+"\t"+"R14" +"\t"+"0"+"\t"+"-",
                   "chr2"+"\t"+"10" +"\t"+"20" +"\t"+"R15" +"\t"+"0"+"\t"+"-",
                   "chr3"+"\t"+"10" +"\t"+"80" +"\t"+"R16" +"\t"+"1"+"\t"+"+",
                   "chr4"+"\t"+"20" +"\t"+"30" +"\t"+"R17" +"\t"+"1"+"\t"+"+",
                   "chr4"+"\t"+"40" +"\t"+"50" +"\t"+"R18" +"\t"+"1"+"\t"+"-",
                   "chr5"+"\t"+"40" +"\t"+"50" +"\t"+"R19" +"\t"+"1"+"\t"+"-"]
    s2_elements = ["chr1"+"\t"+"10" +"\t"+"20" +"\t"+"R21" +"\t"+"0"+"\t"+"+",
                   "chr1"+"\t"+"30" +"\t"+"50" +"\t"+"R22" +"\t"+"1"+"\t"+"-",
                   "chr1"+"\t"+"60" +"\t"+"70" +"\t"+"R23" +"\t"+"0"+"\t"+"+",
                   "chr1"+"\t"+"80" +"\t"+"110"+"\t"+"R24" +"\t"+"4"+"\t"+"+",
                   "chr1"+"\t"+"130"+"\t"+"140"+"\t"+"R25" +"\t"+"0"+"\t"+"+",
                   "chr1"+"\t"+"150"+"\t"+"170"+"\t"+"R26" +"\t"+"1"+"\t"+"-",
                   "chr1"+"\t"+"180"+"\t"+"220"+"\t"+"R27" +"\t"+"0"+"\t"+"-",
                   "chr2"+"\t"+"30" +"\t"+"40" +"\t"+"R28" +"\t"+"0"+"\t"+"-",
                   "chr3"+"\t"+"20" +"\t"+"30" +"\t"+"R29" +"\t"+"0"+"\t"+"-",
                   "chr3"+"\t"+"40" +"\t"+"50" +"\t"+"R210"+"\t"+"0"+"\t"+"-",
                   "chr3"+"\t"+"60" +"\t"+"70" +"\t"+"R211"+"\t"+"0"+"\t"+"+",
                   "chr4"+"\t"+"10" +"\t"+"60" +"\t"+"R212"+"\t"+"0"+"\t"+"+",
                   "chr5"+"\t"+"10" +"\t"+"20" +"\t"+"R213"+"\t"+"1"+"\t"+"-"]
    expect_str =  ["chr1"+"\t"+"40" +"\t"+"50" +"\t"+"X"   +"\t"+"0"+"\t"+"+",
                   "chr1"+"\t"+"60" +"\t"+"70" +"\t"+"X"   +"\t"+"0"+"\t"+"+",
                   "chr1"+"\t"+"80" +"\t"+"90" +"\t"+"X"   +"\t"+"0"+"\t"+"+",
                   "chr1"+"\t"+"100"+"\t"+"110"+"\t"+"X"   +"\t"+"0"+"\t"+"+",
                   "chr1"+"\t"+"160"+"\t"+"170"+"\t"+"X"   +"\t"+"0"+"\t"+"+",
                   "chr1"+"\t"+"180"+"\t"+"190"+"\t"+"X"   +"\t"+"0"+"\t"+"+",
                   "chr1"+"\t"+"200"+"\t"+"210"+"\t"+"X"   +"\t"+"0"+"\t"+"+",
                   "chr3"+"\t"+"20" +"\t"+"30" +"\t"+"X"   +"\t"+"0"+"\t"+"+",
                   "chr3"+"\t"+"40" +"\t"+"50" +"\t"+"X"   +"\t"+"0"+"\t"+"+",
                   "chr3"+"\t"+"60" +"\t"+"70" +"\t"+"X"   +"\t"+"0"+"\t"+"+",
                   "chr4"+"\t"+"20" +"\t"+"30" +"\t"+"X"   +"\t"+"0"+"\t"+"+",
                   "chr4"+"\t"+"40" +"\t"+"50" +"\t"+"X"   +"\t"+"0"+"\t"+"+"]
    input_s1 = [parseBEDString(x) for x in s1_elements]
    input_s2 = [parseBEDString(x) for x in s2_elements]
    expect = [parseBEDString(x) for x in expect_str]
    got = regionsIntersection(input_s1, input_s2)
    if debug :
      sys.stderr.write("expect:\n")
      sys.stderr.write("\n".join([str(x) for x in expect]) + "\n")
      sys.stderr.write("got:\n")
      sys.stderr.write("\n".join([str(x) for x in got]) + "\n")
    self.assertEqual(expect, got)
      
  def testSizeOfOverlap(self):
    ## region 2 entirely inside region 1 -- ans = size of region 2
    region1 = GenomicInterval("chr1", 10, 20, "read", 1, "f")
    region2 = GenomicInterval("chr1", 10, 20, "read", 1, "f")
    self.assertEqual(region1.sizeOfOverlap(region2), len(region2))
    self.assertEqual(region1.sizeOfOverlap(region2), 
                     region2.sizeOfOverlap(region1))
  
  def testSubtract(self):
    debug = False
    a = GenomicInterval("chr1",10,100)
    b1 = GenomicInterval("chr1",1,8)
    b2 = GenomicInterval("chr1",15,20)
    b3 = GenomicInterval("chr1",50,60)
    b4 = GenomicInterval("chr1",90,120)
    c1 = GenomicInterval("chr2",30,40)
    
    res = a.subtract([b1,b2,b3,b4,c1])
    expect = [GenomicInterval("chr1",10,15),
              GenomicInterval("chr1",20,50),
              GenomicInterval("chr1",60,90)]
    res.sort()
    expect.sort()
    if debug :
      sys.stderr.write("\nDebugging BED subtraction test\n")
      sys.stderr.write("expect\n" + "\n".join([str(x) for x in expect]) + "\n")
      sys.stderr.write("got\n" + "\n".join([str(x) for x in res]) + "\n")
    assert(res == expect)
  
  def testDistance(self):
    """
      @summary: test that calculating the distance between two GenomicInterval
                objects succeeds when they're on the same chromosome and gives
                the right answer, but fails when they're on different 
                chromosomes
      @todo: needs to be implemented
    """
    pass
  
  def testSignedDistance(self):
    """
      @summary: test that calculating the signed distance between two 
                GenomicInterval objects succeeds when they're on the same 
                chromosome and gives the right answer, but fails when they're  
                on different chromosomes
      @todo: needs to be implemented
    """
    pass
  
  def testIntersects(self):
    """
      @summary: test that intervals on different chroms don't intersect even
                when they have the same co-ordinates. Test that intervals 
                on the same chrom don't intersect when only the end co-ordinate
                overlaps the start of the other. 
      @todo: test needs to be implemented
    """
    pass


if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])

