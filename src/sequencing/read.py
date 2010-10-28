#!/usr/bin/python

""" 
  Date of Creation: 20th May 2010     
                       
  Description:   defines fasta and fastq read classes 

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
  History:       2nd Sep 2010 -- Philip Uren
                   * added this header
                   * added trim methods
                 4th Sep 2010 -- Philip Uren
                   * added effective length method 
                   * added simple unit tests for length methods
                 13th Sep 2010 -- Philip Uren
                   * moved common code for fastq/a to superclass 
                   * removed inline c code (for now..)
                 17th Sep 2010 -- Philip Uren
                   * added formattedString method
                   * added unit test for formattedString method 
                 29th Sep 2010 -- Philip Uren
                   * fixed bug in trimRight 
  
  TODO:          Move these classes into the associate fastq and fasta files? 
                 truncate method should move up class hierarchy -- problem with dynamic linking? 
                 unit tests missing  
"""

import unittest, sys
from fastread import FastRead

class FastqReadError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)
  
class FastaRead(FastRead):
  def __init__(self, seqName, seqData = "", useMutableString = False):
    FastRead.__init__(self, seqName, seqData, useMutableString)
  def __str__(self):
    """
      DESCPT: return string representation of the read
    """
    return ">" + self.sequenceName + "\n" + str(self.sequenceData) 
  
  def formattedString(self, width):
    """
      DESCPT: get formatted version of the seq where no 
              row of sequence data exceeds a given length
    """
    res = ">" + self.sequenceName + "\n"
    for i in range(0,len(self.sequenceData), width) :
      res += self.sequenceData[i:i+width]
      if i + width < len(self.sequenceData) : res += "\n" 
    return res
        
class FastqRead(FastRead):
  def __init__(self, seqName, seqData=None, 
               seqQual=None, useMutableString=False):
    FastRead.__init__(self, seqName, seqData, useMutableString)
    self.sequenceQual = seqQual
      
    # for quality scores 
    self.LOWSET_SCORE = 64
    self.HIGHEST_SCORE = 104
        
  def __eq__(self, read):
    FastRead.__eq__(read) and self.sequenceQual == read.sequenceQual
            
  def __ne__(self, read):
    FastRead.__ne__(read) and self.sequenceQual != read.sequenceQual
    
  def truncate(self, size):
    self.trimRight(len(self) - size)
        
  def trimRight(self, amount):
    """
      @summary: trim the read by removing <amount> nucleotides
                from the 3' end (right end)
    """
    self.sequenceData = self.sequenceData[:-amount]
    self.sequenceQual = self.sequenceQual[:-amount]
      
  def trimLeft(self, amount):
    self.sequenceData = self.sequenceData[amount:]
    self.sequenceQual = self.sequenceQual[amount:]
    
  def getRelativeQualityScore(self, i):
    val = self.sequenceQual[i]
    return (ord(val) - self.LOWSET_SCORE) / float (self.HIGHEST_SCORE - self.LOWSET_SCORE)

  def split(self, point = None):
    """ returns two Read objects which correspond to the split of this read """
    if point == None :
      point = len(self)/2
    
    r1 = FastqRead(self.sequenceName + ".1", 
              self.sequenceData[:point], 
              self.sequenceQual[:point])
    r2 = FastqRead(self.sequenceName + ".2", 
              self.sequenceData[point:],
              self.sequenceQual[point:])
    return r1,r2
      
  def merge(self, other, forceMerge = False):
    if self.sequenceName != other.sequenceName and not forceMerge :
      raise FastqReadError("cannot merge " + self.sequenceName + " with " + other.sequenceName + " -- different sequence names")
    
    name = self.sequenceName
    seq = self.sequenceData + other.sequenceData
    qual = self.sequenceQual + other.sequenceQual
    
    return FastqRead(name, seq, qual)
                  
  def __str__(self):
    return "@" + self.sequenceName + "\n" + self.sequenceData +\
           "\n" + "+" + self.sequenceName + "\n" + self.sequenceQual 
    

class ReadUnitTests(unittest.TestCase):
  """
    unit tests for read.py 
  """
    
  def SetUp(self):
    pass
  
  def testFormattedString(self):
    r = FastaRead("name", "ATCGATCGATCGATCTCGA")
    expect = ">name\n" +\
             "ATCGA\n" +\
             "TCGAT\n" +\
             "CGATC\n" +\
             "TCGA"
    got = r.formattedString(width=5)
    self.assertTrue(got == expect)
    
    # make sure this also works with a mutable underlying sequence
    r = FastaRead("name", "ATCGATCGATCGATCTCGA", useMutableString = True)
    expect = ">name\n" +\
             "ATCGA\n" +\
             "TCGAT\n" +\
             "CGATC\n" +\
             "TCGA"
    got = r.formattedString(width=5)
    self.assertTrue(got == expect)
    
  
if __name__ == "__main__":
    unittest.main()