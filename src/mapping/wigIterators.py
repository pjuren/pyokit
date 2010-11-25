#!/usr/bin/python

""" 
  Date of Creation: 11th November 2010     
                       
  Description:   Iterators for processing wiggle format data streams

  Copyright (C) 2010  
  University of Southern California,
  Philip J. Uren,
  Andrew D. Smith
  
  Authors: Philip J. Uren, Andrew D. Smith
  
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
  History:       16th November 2010 -- Philip Uren
                   * added sorting parameter to wigIterator
  
  TODO:          * Class and method headers are missing
"""

ITERATOR_SORTED_START = 1

import sys, unittest
from util.fileUtils import openFD, getFDName
from mapping.wig import wigElementFromString, WigError
from testing.dummyfiles import DummyInputStream, DummyOutputStream
from util.fileUtils import linesInFile
from util.progressIndicator import ProgressIndicator

def wigIterator(fd, verbose = False, sortedby = None):
  """
    @param sortedBy: if not None, should be one of ITERATOR_SORTED_BY_START
                     indicating an order that the input stream must be 
                     sorted in
    @raise WigError: if sortedBy is set and stream is not sorted
  """
  if verbose :
    try :
      totalLines = linesInFile(fd)
      pind = ProgressIndicator(totalToDo = totalLines, 
                               messagePrefix = "completed", 
                               messageSuffix = "of processing " +\
                                               getFDName(fd))
    except AttributeError :
      sys.stderr.write("WigIterator -- warning: " +\
                       "unable to show progress for stream")
      verbose = False   
      
  chromsSeen = set()
  prev = None
  
  fh = openFD(fd)
  for line in fh : 
    if verbose :
      pind.done += 1
      pind.showProgress()
    
    line = line.strip()
    if line == "" : continue 
    e = wigElementFromString(line)
    
    # on same chrom as the prev item, make sure order is right
    if prev != None and sortedby != None and e.chrom == prev.chrom :
      if sortedby == ITERATOR_SORTED_START and prev.start > e.start :
        raise BEDError("bed file " + filehandle.name +\
                       " not sorted by start index - saw item " +\
                       str(prev) + " before " + str(e))
    
    # starting a new chrom.. make sure we haven't already seen it
    if prev != None and prev.chrom != e.chrom :
      if (sortedby == ITERATOR_SORTED_START) and\
         (e.chrom in chromsSeen or prev.chrom > e.chrom) :
        raise BEDError("BED file " + filehandle.name +\
                       " not sorted by chrom")
      chromsSeen.add(e.chrom) 
    
    # all good..
    yield e
    prev = e
    
    
class WigIteratorUnitTests(unittest.TestCase):
  """
    Unit tests for wigIterator 
  """
  
  def setUp(self):
    pass

  def testWigIterator(self):
    debug = False
    wigIn = "chr1 \t 1 \t 2 \t 5\n" +\
            "chr1 \t 40 \t 50 \t 3\n" +\
            "chr2 \t 30 \t 34 \t 2\n" +\
            "chr4 \t 30 \t 60 \t 6\n"
    expect = ["chr1\t1\t2\t5",
              "chr1\t40\t50\t3",
              "chr2\t30\t34\t2",
              "chr4\t30\t60\t6"] 
    
    out = []
    for element in wigIterator(DummyInputStream(wigIn)) :
      out.append(str(element))
    
    out.sort()
    expect.sort()     

    if debug :
      print "out ------ "
      for l in out : print l
      print "-------"
      print "expect ------ "
      for l in expect : print l
      print "-------"
    
    self.assertTrue(out == expect)
    

if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])