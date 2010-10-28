#!/usr/bin/python

""" 
  Date of Creation: 13th September 2010     
                       
  Description:   Defines superclass for fastq and fasta reads

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
  History:       16th September 2010 -- Philip Uren
                   * added mask regions methods
                   * added unit tests for masking methods
                   * added copy constructor, and unit tests
                   * fixed bug in similarity 
                   * fixed bug in clipThreePrime
                 17th September -- Philip Uren
                   * added mutable string class
                   * allow underlying string rep for sequence to be 
                     mutable string
                   * improved maskRegion efficiency when using 
                     mutable strings for sequence  
                   * fixed behaviour of mutable string for slices
                 27th October -- Philip Uren 
                   * moved code to smithlab_py
                   * fixed import statements to reflect change in
                     project structure 
  
  
  TODO:          Expand unit tests
"""

import unittest, sys, os
from util.progressIndicator import ProgressIndicator 

class MutableString :
  def __init__(self, strng):
    self.list = []
    for ch in strng :
      self.list.append(ch)
  
  def __getitem__(self, indx):
    if isinstance(indx, slice):
      return "".join(self.list.__getitem__(indx))
    return self.list[indx]
  
  def __setitem__(self, k, v):
    self.list[k] = v
  
  def __len__(self):
    return len(self.list)
  
  def __eq__(self, other):
    if type(other).__name__ == "str" :
      other = [x for x in other]
    else :
      try :
        other = other.list
      except : return false 
    return self.list == other
  
  def __ne__(self, other):
    if type(other).__name__ == "str" :
      other = [x for x in other]
    else :
      try :
        other = other.list
      except : return false 
    return self.list != other
  
  def __str__(self):
    return "".join(self.list)
  

class FastRead:
  def __init__(self, seqName, seqData, useMutableString = False):
    """
      DESCPT: Constructor 
    """
    self.sequenceName = seqName
    if useMutableString : 
      self.sequenceData = MutableString(seqData)
    else :
      self.sequenceData = seqData
    self.mutableString = useMutableString
    
  def copy(self):
    """
      DESCPT: Copy constructor 
    """
    return FastRead(self.sequenceName, self.sequenceData, self.mutableString)
    
  def percentNuc(self, nuc):
    """
      DESCPT: return the percentage of the sequence which is equal 
              to the passed nucleotide
      PARAMS: nuc   -- count occurences of this nuc  
      RETURN: float -- percentage
    """
    count = reduce(lambda x,y: x+1 if y==nuc else x, self.sequenceData, 0) 
    return count / float(len(self.sequenceData))
  
  def similarity(self, self_start, self_end, other_start, other_end, other):
    """
      DESCPT: Count of the number of places where this[start,end] is equal to
              other[o_start, o_end]
    """ 
    assert(self_end - self_start == other_end - other_start)
    count = 0
    for i in range(0, self_end - self_start + 1) :
      if self.sequenceData[self_start + i] == other.sequenceData[other_start + i] : 
        count += 1
    return count
    
  def __len__(self):
    """
      DESCPT: length of the read, defined as the length of its sequence data
    """
    return len(self.sequenceData)
  
  def effectiveLength(self):
    """ 
      DESCPT: disregarding N's, how long is this sequence? 
    """
    return len([nuc for nuc in self.sequenceData 
                    if nuc != "N" and nuc != "n"])
  
  def __eq__(self, read):
    """
      DESCRP: return true if this read is equal to passed parameter
    """
    if read == None : return False
    return  self.sequenceData == read.sequenceData and\
            self.sequenceName == read.sequenceName
              
  def __ne__(self, read):
    """
      DECPRP: return true if this read is not equal to the passed parameter  
    """
    if read == None : return True
    return  self.sequenceData != read.sequenceData or\
            self.sequenceName != read.sequenceName
            
  def nsLeft(self, amount):
    """
      DESCPT: replace leftmost <amount> bases by Ns
    """
    self.sequenceData = (amount * "N") + self.sequenceData[amount:]
    
  def nsRight(self, amount):
    """
      DESCPT: replace rightmost <amount> bases by Ns
    """
    self.sequenceData = self.sequenceData[:-amount] + (amount * "N")
    
  def maskRegion(self, region):
    """
      DESCPT: Replace nucleotides in this read in the regions 
              given by Ns 
      PARAMS: region -- any object with .start and .end attributes
                        co-ords are zero based and inclusive of both 
                        end points 
      RAISES: FastreadError -- if region specifies nucleotides not present in
                               this read
    """
    if region.start < 0 or region.end < 0 or \
       region.start > len(self) or region.end > len(self) :
      raise FastreadError("cannot mask region " + str(region.start) + " to " +\
                          str(region.end) + " in " + self.sequenceName + ". " +\
                          "Region specifies nucleotides not present in " +\
                          "this read. Valid range would have been 0 -- " +\
                          str(len(self)))
      
    if self.mutableString : 
      for i in range(region.start, region.end + 1) : self.sequenceData[i] = 'N'
    else :
      self.sequenceData = "".join([self.sequenceData[:region.start],
                          ("N" * (region.end - region.start + 1)),
                          self.sequenceData[region.end+1:]])
                        
  def maskRegions(self, regions, verbose = False):
    """
      DESCPT: Mask the given regions in this read
      PARAMS: verbose -- print status messages to stderr if True
    """
    if verbose:
      pind = ProgressIndicator(totalToDo = len(regions), 
                               messagePrefix = "completed", 
                               messageSuffix = "of masking regions in " +\
                                                self.sequenceName)
    for region in regions :
      self.maskRegion(region)
      if verbose : 
        pind.done += 1
        pind.showProgress()
      
      
  def split(self, point = None):
    """ 
      DESCPT: Split this read into two halves
      PARAMS: point -- defines the split point, if None then the centre is used
      RETURN: two FastRead objects -- one for each side  
    """
    if point == None :
      point = len(self)/2
    
    r1 = FastqRead(self.sequenceName + ".1", 
                   self.sequenceData[:point])
    r2 = FastqRead(self.sequenceName + ".2", 
                   self.sequenceData[point:])
    
    return r1,r2
  
  def truncate(self, newLength):
    """
      DESCRP: truncate the read so it is only <newLength> nucleotides long
    """
    return FastRead(self.sequenceName, self.sequenceData[:10])
  
  def clipThreePrime(self, seq, mm_score):
    """
      DESCPT: Clip a sequence from the 3' end of the read -- we assume sequence to
              be clipped will always begin somewhere in this sequence, but may 
              not be fully contained. If found, replaced with Ns
      PARAMS: seq -- sequence to be clipped
              mm_score -- the number of matching bases needed to consider a hit,
                          mm_score = len(adaptor) would be 100% match 
    """
    lim = mm_score - 1
    other_end = len(seq) - 1
    other_start = 0
  
    for i in range(len(self.sequenceData) - 1, lim - 1, -1) :
      self_end = i
      self_start = i - (len(seq) - 1)
      
      if self_start < 0 : 
        self_start = 0
        other_start = other_end - self_end
      
      score = self.similarity(self_start, self_end, other_start, other_end, seq)
      if (score >= mm_score) :
        self.nsRight(len(seq) + (len(self) - i) - 1)
        break
  
  def clipAdaptor(self, adaptor):
    """
      DESCPT: Clip adaptor sequence from this read. We assume it's in the 
              3' end
      PARAMS: adaptor -- sequence to look for. We only use first 10 bases
                         must be a full FastRead object, not just string
    """
    missmatches = 2
    adaptor = adaptor.truncate(10)
    self.clipThreePrime(adaptor, len(adaptor) - missmatches)
    
  def containsAdaptor(self, adaptor):
    """
      DESCPT: Does this sequence contain adaptor contamination? 
              We assume adaptor is in 3' end
      PARAMS: adaptor -- sequence to look for. must be a full FastRead 
                         object, not just string
      RETURN: bool -- true if there is an occurence of <adaptor>, 
              false otherwise
    """
    origSeq = self.sequenceData
    self.clipAdaptor(adaptor)
    res = False
    if self.sequenceData != origSeq : res = True
    self.sequenceData = origSeq
    return res
  
  def isPolyA(self):
    """
      DESCRP: Is this sequence a polyA? based on having > 90% A in read
    """
    return self.percentNuc("A") >= 0.9
  
  def isPolyT(self):
    """
      DESCRP: Is this sequence a polyT? based on having > 90% A in read
    """
    return self.percentNuc("T") >= 0.9
  
  def isLowQuality(self):
    """
      DESCRP: Is this sequence low quality? based on having > 10% N in read
    """
    return self.percentNuc("N") >= 0.1
  
  def maskMatch(self, mask):
    """
      DESCPT: Determine whether this sequence matches the given mask.
              Ns in the mask are considered to match anything in the
              sequence -- all other chars must match exactly 
      PARAMS: mask -- string to match against
      RETURN: True if the mask matches at all places, otherwise false
    """
    if len(mask) > len(self.sequenceData) : return False
    lim = len(mask)
    for i in range(0,lim):
      if mask[i] == "N" or mask[i] == "n" : continue
      if mask[i] != self.sequenceData[i] : return False
    return True
  
  
class FastReadUnitTests(unittest.TestCase):
  """
    Unit tests for fastread 
  """
    
  def testClipadaptor(self):
    pass
    input =   FastRead("name", "ACTGCTAGCGATCGACT")
    adaptor = FastRead("adap",       "AGCGATAGACT")
    expect =  FastRead("name", "ACTGCTNNNNNNNNNNN")
    input.clipAdaptor(adaptor)
    got = input
    self.assertTrue(expect == got)
    
  def testNsLeft(self):
    input =   FastRead("name", "ACTGCTAGCGATCGACT")
    expect =  FastRead("name", "NNNNNTAGCGATCGACT")
    input.nsLeft(5)
    got = input
    self.assertTrue(expect == got)
    
  def testNsRight(self):
    input =   FastRead("name", "ACTGCTAGCGATCGACT")
    expect =  FastRead("name", "ACTGCTAGCGATNNNNN")
    input.nsRight(5)
    got = input
    self.assertTrue(expect == got)
    
  def testLengths(self):
    input =   FastRead("name", "ACTNCTANCGATNNACT")
    self.assertTrue(len(input) == 17)
    self.assertTrue(input.effectiveLength() == 13)
    
  def testMaskRegion(self):
    class TestRegion:
      def __init__(self, s, e):
        self.start = s
        self.end = e
        
    input =   FastRead("name", "ACTNCTANCGATNNACT")
    expect =  FastRead("name", "ANNNNTANCGATNNACT")
    out = input.copy()
    out.maskRegion(TestRegion(1,4))
    self.assertTrue(expect == out)
    
  def testCopyConstructor(self):
    one = FastRead("name", "ACTNCTANCGATNNACT")
    two = one.copy()
    self.assertTrue(one == two)
    one.sequenceName = "old"
    self.assertTrue(one != two)
    
  
if __name__ == "__main__":
    unittest.main()
