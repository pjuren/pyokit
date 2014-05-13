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
                 28th October -- Philip Uren
                   * moved fastq read classes to fastqread module 
                   * fixed bugs in eq and ne --> added unit tests for eq and ne
  
  TODO:           
                 * truncate method should move up class hierarchy -- problem with 
                   dynamic linking? 
                 * unit tests missing  
                 * method headers missing 
"""

import unittest, sys
from fastread import FastRead

# this is the offset needed to shift Sanger fastq quality scores into solexa 
SANGER_SOLEXA_OFFSET = 31

class FastqReadError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)
        
class FastqRead(FastRead):
  def __init__(self, seqName, seqData=None, 
               seqQual=None, useMutableString=False):
    FastRead.__init__(self, seqName, seqData, useMutableString)
    self.sequenceQual = seqQual
      
    # for quality scores 
    self.LOWSET_SCORE = 64
    self.HIGHEST_SCORE = 104
        
  def __eq__(self, read):
    return FastRead.__eq__(self, read) and \
           self.sequenceQual == read.sequenceQual
            
  def __ne__(self, read):
    return FastRead.__ne__(self, read) or \
           self.sequenceQual != read.sequenceQual
    
  def truncate(self, size):
    self.trimRight(len(self) - size)
    
  def qualityToSolexa(self):
    """
      @summary: convert quality data from sanger to solexa format. Note that
                no checking is done to make sure the data was originally in 
                sanger format; if it wasn't, the result will be junk
    """
    newqual = ""
    for val in self.sequenceQual :
      newqual += chr(ord(val) + SANGER_SOLEXA_OFFSET)
    self.sequenceQual = newqual 
        
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
    return (ord(val) - self.LOWSET_SCORE) / float (self.HIGHEST_SCORE - 
                                                   self.LOWSET_SCORE)

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
      raise FastqReadError("cannot merge " + self.sequenceName + " with " +\
                           other.sequenceName + " -- different sequence names")
    
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
  
  def testeq(self):
    r1 = FastqRead("s1","ACTGCT","BBBBBB")
    r2 = FastqRead("s1","ACTGCT","BBBBBB")
    r3 = FastqRead("s1","ACTGCT","BBBfBB")
    r4 = FastqRead("s1","ACCGCT","BBBBBB")
    r5 = FastqRead("s2","TCTGCT","fBBBBB")
    r6 = FastqRead("s3","CCCCCC","fBBBBB")
    r7 = FastqRead("s1","CCCCCC","fBBfBB")
    r7 = FastqRead("s6","CCCCCC","fBBfBB")
    
    self.assertTrue((r1 == r2) == True)  # same name, same seq, same qual
    self.assertTrue((r1 == r3) == False) # same name, same seq, diff qual
    self.assertTrue((r1 == r4) == False) # same name, diff seq, same qual
    self.assertTrue((r3 == r4) == False) # same name, diff seq, diff qual
    self.assertTrue((r1 == r5) == False) # diff name, diff seq, diff qual
    self.assertTrue((r5 == r6) == False) # diff name, diff seq, same qual
    self.assertTrue((r6 == r7) == False) # diff name, same seq, diff qual
    self.assertTrue((r3 == r4) == False) # diff name, same seq, same qual
    
    self.assertTrue((r1 != r2) == False)  # same name, same seq, same qual
    self.assertTrue((r1 != r3) == True) # same name, same seq, diff qual
    self.assertTrue((r1 != r4) == True) # same name, diff seq, same qual
    self.assertTrue((r3 != r4) == True) # same name, diff seq, diff qual
    self.assertTrue((r1 != r5) == True) # diff name, diff seq, diff qual
    self.assertTrue((r5 != r6) == True) # diff name, diff seq, same qual
    self.assertTrue((r6 != r7) == True) # diff name, same seq, diff qual
    self.assertTrue((r3 != r4) == True) # diff name, same seq, same qual

  
if __name__ == "__main__":
    unittest.main()