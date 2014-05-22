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
                 28th October 2010 -- Philip Uren
                   * Moved fasta classes into fastaread module 
  
  TODO:           
                 * unit tests missing  
"""

import unittest, sys
from pyokit.sequencing.fastread import FastRead

class FastaReadError(Exception):
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