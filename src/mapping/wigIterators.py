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
                 23rd December 2010 -- Philip Uren
                   * added pairedWigIterator
                 18th October 2011 -- Philip Uren
                   * added fixed step wig file iterator 
                   * allowed wigIterator to automatically determine file type
  
  TODO:          * Class, method and function comment headers are incomplete 
"""

ITERATOR_SORTED_START = 1

import sys, unittest, os
from util.fileUtils import openFD, getFDName
from mapping.wig import wigElementFromString, WigError, WigElement
from testing.dummyfiles import DummyInputStream, DummyOutputStream
from util.fileUtils import linesInFile
from util.progressIndicator import ProgressIndicator

def wigIterator(fd, verbose=False, sortedby=None):
  # peak at the first line to see if it's a regular wig, or 
  # fixed-step wig
  fh = openFD(fd)
  at = fh.tell()
  line = None
  while(line==None) :
    l = fh.readline()
    if l.strip() != "" : line = l
  
  fh.seek(at) 
  if line.split()[0] == "fixedStep" : 
    return fixedWigIterator(fd,verbose,sortedby)
  else :
    return regularWigIterator(fd,verbose,sortedby)
  

def regularWigIterator(fd, verbose = False, sortedby = None):
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
    try :
      e = wigElementFromString(line)
    except :
      raise WigError("failed parsing '" + line + "' as wig element")
    
    # on same chrom as the prev item, make sure order is right
    if prev != None and sortedby != None and e.chrom == prev.chrom :
      if sortedby == ITERATOR_SORTED_START and prev.start > e.start :
        raise WigError("bed file " + fd.name +\
                       " not sorted by start index - saw item " +\
                       str(prev) + " before " + str(e))
    
    # starting a new chrom.. make sure we haven't already seen it
    if prev != None and prev.chrom != e.chrom :
      if (sortedby == ITERATOR_SORTED_START) and\
         (e.chrom in chromsSeen or prev.chrom > e.chrom) :
        raise WigError("BED file " + fd.name +\
                       " not sorted by chrom")
      chromsSeen.add(e.chrom) 
    
    # all good..
    yield e
    prev = e
    
def fixedWigIterator(fd, verbose=False, sortedby = None):
  """
    @summary: 
  """
  fh = openFD(fd)
  if verbose :
    try :
      pind = ProgressIndicator(totalToDo = os.path.getsize(fh.name), 
                                     messagePrefix = "completed", 
                                     messageSuffix = "of processing " +\
                                                      fh.name)
    except AttributeError :
      sys.stderr.write("WigIterator -- warning: " +\
                       "unable to show progress for stream")
      verbose = False   
      
  chromsSeen = set()
  prev = None
  
  #NUMBERS = set(['1','2','3','4','5','6','7','8','9','0','.'])
  currentChrom, at, step = None, None, None
  for line in fh : 
    line = line.strip()
    if line == "" : continue
    
    if line[0] == 't' or line[0] == 'f' :
      parts = line.split()
      if parts[0] == "track" : continue 
      elif parts[0] == "fixedStep" :
        currentChrom = parts[1].split("=")[1]
        at = int(parts[2].split("=")[1])
        step = int(parts[3].split("=")[1])
    else :
      val = float(line)
      e = WigElement(currentChrom, at, at+step, val)
    
      # on same chrom as the prev item, make sure order is right
      if prev != None and sortedby != None and e.chrom == prev.chrom :
        if sortedby == ITERATOR_SORTED_START and prev.start > e.start :
          raise WigError("bed file " + fd.name +\
                         " not sorted by start index - saw item " +\
                         str(prev) + " before " + str(e))
    
      # starting a new chrom.. make sure we haven't already seen it
      if prev != None and prev.chrom != e.chrom :
        if (sortedby == ITERATOR_SORTED_START) and\
           (e.chrom in chromsSeen or prev.chrom > e.chrom) :
          raise WigError("BED file " + fd.name +\
                         " not sorted by chrom")
        chromsSeen.add(e.chrom) 
    
      # all good..
      yield e
      prev = e
      at += step
      if verbose :
        pind.done = fh.tell()
        pind.showProgress()

    
def pairedWigIterator(wgs1, wgs2, missingVal = 0, verbose = False, debug = False):
  """
    @summary: iterate over two wig streams, returning matching pairs from each
              file. If a match is missing for a given item in one file then 
              it is assumed that item has the value <missingVal>, which by
              default is 0. 
    @note: streams must be sorted by chrom and start coordinate.
    @note: size of each wig element must be the same as its pair, else an 
           exception is raised
  """
  BOTH, FIRST, SECOND = 0, 1, 2
  wigi1 = wigIterator(wgs1, verbose = verbose, sortedby = ITERATOR_SORTED_START)
  wigi2 = wigIterator(wgs2, verbose = verbose, sortedby = ITERATOR_SORTED_START)
  citem1, citem2 = None, None
  pitem1, pitem2 = None, None
  next = BOTH
  
  while True :
    # try to get items
    if next == BOTH or next == FIRST :
      try :
        citem1 = wigi1.next()
      except StopIteration:
        citem1 = None
    if next == BOTH or next == SECOND : 
      try :
        citem2 = wigi2.next()
      except StopIteration:
        citem2 = None
    
    if debug : print "got items: " + str(citem1) + " and " + str(citem2)
    
    # if both streams are exhausted, we're done
    if citem1 == None and citem2 == None : break
    
    # if we got two things that match, yield them and remember that next time 
    # we want something from both streams
    if citem1 != None and citem2 != None and \
       citem1.start == citem2.start and citem1.end == citem2.end and \
       citem1.chrom == citem2.chrom :
      next = BOTH
      if debug : print "output: " + str(citem1) + " and " + str(citem2) + "\n"
      yield citem1, citem2
      continue 
    
    # we got non-matching items, the smaller one can be yielded with a dummy
    # wig element for the missing value, the larger needs to be kept so we 
    # can check for it's match further along the stream
    if citem2 == None or (citem1 != None and citem1.before(citem2)) :
      next = FIRST
      p1 = citem1 
      p2 = WigElement(citem1.chrom, citem1.start, citem1.end, missingVal)
    elif citem1 == None or (citem2 != None and citem2.before(citem1)) :
      next = SECOND
      p1 = WigElement(citem2.chrom, citem2.start, citem2.end, missingVal)
      p2 = citem2
    if debug : print "output: " + str(p1) + " and " + str(p2) + "\n"
    yield p1, p2
    
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
    
  def testFixedWigIterator(self):
    debug = False
    wigIn = "fixedStep chrom=chr1 start=1 step=1\n" +\
            "5\n" +\
            "3\n" +\
            "fixedStep chrom=chr2 start=30 step=2\n" +\
            "2\n" +\
            "4\n" +\
            "fixedStep chrom=chr4 start=10 step=1\n" +\
            "6\n" +\
            "0.5\n"
    expect = ["chr1" + "\t" +  "1" + "\t" +  "2" +"\t" + "5",
              "chr1" + "\t" +  "2" + "\t" +  "3" +"\t" + "3",
              "chr2" + "\t" + "30" + "\t" + "32" +"\t" + "2",
              "chr2" + "\t" + "32" + "\t" + "34" +"\t" + "4",
              "chr4" + "\t" + "10" + "\t" + "11" +"\t" + "6",
              "chr4" + "\t" + "11" + "\t" + "12" +"\t" + "0.5"] 
    
    out = []
    for element in fixedWigIterator(DummyInputStream(wigIn)) :
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
    
  def testPairedWigIterator(self):
    debug = False
    wigIn1 = "chr1 \t 10 \t 20 \t 01\n" +\
             "chr1 \t 20 \t 30 \t 02\n" +\
             "chr1 \t 30 \t 40 \t 03\n" +\
             "chr2 \t 40 \t 50 \t 04\n" +\
             "chr3 \t 70 \t 80 \t 05\n" +\
             "chr4 \t 10 \t 20 \t 06\n" 
    wigIn2 = "chr1 \t 10 \t 20 \t 99\n" +\
             "chr1 \t 30 \t 40 \t 98\n" +\
             "chr2 \t 40 \t 50 \t 97\n" +\
             "chr2 \t 50 \t 60 \t 96\n" +\
             "chr3 \t 70 \t 80 \t 95\n" +\
             "chr5 \t 10 \t 20 \t 94\n" 
    expect = [("chr1\t10\t20\t1", "chr1\t10\t20\t99"),
              ("chr1\t20\t30\t2", "chr1\t20\t30\t-1"),
              ("chr1\t30\t40\t3", "chr1\t30\t40\t98"),
              ("chr2\t40\t50\t4", "chr2\t40\t50\t97"),
              ("chr2\t50\t60\t-1", "chr2\t50\t60\t96"),
              ("chr3\t70\t80\t5", "chr3\t70\t80\t95"),
              ("chr4\t10\t20\t6", "chr4\t10\t20\t-1"),
              ("chr5\t10\t20\t-1", "chr5\t10\t20\t94")]
    
    out = []
    for e1,e2 in pairedWigIterator(DummyInputStream(wigIn1), 
                                   DummyInputStream(wigIn2), 
                                   missingVal = -1, debug = debug) :
      out.append((str(e1), str(e2)))
    
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