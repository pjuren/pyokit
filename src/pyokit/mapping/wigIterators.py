#!/usr/bin/python

"""
  Date of Creation: 11th November 2010
  Description:      Iterators for processing wiggle format data streams

  Copyright (C) 2010-2014
  Philip J. Uren

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

ITERATOR_SORTED_START = 1
ITERATOR_SORTED_END = 2

# standard library imports
import sys, unittest, os

# Pyokit imports
from pyokit.datastruct.genomicInterval import GenomicInterval, parseWigString
from pyokit.util.fileUtils import openFD, getFDName
from pyokit.testing.dummyfiles import DummyInputStream, DummyOutputStream
from pyokit.util.fileUtils import linesInFile
from pyokit.util.progressIndicator import ProgressIndicator
from operator import itemgetter, attrgetter

def wigIterator(fd, verbose=False, sortedby=None, scoreType=int):
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
    return regularWigIterator(fd,verbose,sortedby,scoreType=scoreType)


def regularWigIterator(fd, verbose = False, sortedby = None, scoreType=int):
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
    e = parseWigString(line, scoreType=scoreType)

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

def fixedWigIterator(fd, verbose=False, sortedby = None, scoreType = int):
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
      e = GenomicInterval(currentChrom, at, at+step, None,
                          val, scoreType = scoreType)

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


def pairedWigIterator(inputStreams, mirror=False, mirrorScore=None,
                      ignoreScore=True, sortedby=ITERATOR_SORTED_END,
                      scoreType=int, verbose = False, debug = False):
  """
    @summary: iterate over multiple wig streams, and yield a list of wig
              elements that match for each location (locations with 0 matching
              items are skipped)
    @param inputStrams: TODO
    @param mirror:      TODO
    @param mirrorScore: TODO
    @param ignoreScore: Don't consider score when determining if two elements
                        are equal
    @param sortedby:    TODO
    @param scoreType:   TODO
    @param verbose:     TODO
    @param debug:       TODO
    @note: streams must be sorted -- sortedby parameter determines what sorting
           order will be acceptable
  """
  # let's build our sorting order...
  sortOrder = ["chrom"]
  if sortedby == ITERATOR_SORTED_START :
    sortOrder.append("start")
    sortOrder.append("end")
  elif sortedby == ITERATOR_SORTED_END :
    sortOrder.append("end")
    sortOrder.append("start")
  if not ignoreScore : sortOrder.append("score")
  keyFunc = attrgetter(*sortOrder)

  def next(iterator):
    """ little internal function to return the next item, or None """
    try : return iterator.next()
    except StopIteration : return None

  wIterators = [wigIterator(fh, verbose=verbose, sortedby=sortedby,
                            scoreType=scoreType) for fh in inputStreams]
  elements = [next(it) for it in wIterators]

  while True :
    assert(len(elements) >= 2)
    if None not in elements and len(set([keyFunc(x) for x in elements])) == 1 :
      # All equal -- yield and move on for all streams
      yield [e for e in elements]
      elements = [next(it) for it in wIterators]
    else :
      # something wasn't equal.... find the smallest thing, it's about
      # to drop out of range...
      minElement = min([x for x in elements if x != None], key=keyFunc)
      minIndices = [i for i in range(0, len(elements))
                    if elements[i] != None and
                       keyFunc(elements[i]) == keyFunc(minElement)]
      if mirror :
        # mirror the min item for any streams in which it doesn't match
        score = minElement.score if mirrorScore == None else mirrorScore
        yield [elements[i] if i in minIndices
               else GenomicInterval(minElement.chrom, minElement.start,
                                    minElement.end, None, score)
               for i in range(0, len(elements))]

      # move the smallest element onwards now, we're done with it
      for index in minIndices : elements[index] = next(wIterators[index])

    # stop once all strams are exhuasted
    if reduce(lambda x,y:x and y, [e == None for e in elements]) : break


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
    expect = ["chr1" + "\t" +  "1" + "\t" +  "2" +"\t" + "5.0",
              "chr1" + "\t" +  "2" + "\t" +  "3" +"\t" + "3.0",
              "chr2" + "\t" + "30" + "\t" + "32" +"\t" + "2.0",
              "chr2" + "\t" + "32" + "\t" + "34" +"\t" + "4.0",
              "chr4" + "\t" + "10" + "\t" + "11" +"\t" + "6.0",
              "chr4" + "\t" + "11" + "\t" + "12" +"\t" + "0.5"]

    out = []
    for element in fixedWigIterator(DummyInputStream(wigIn), scoreType=float) :
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
    for e1,e2 in pairedWigIterator([DummyInputStream(wigIn1),
                                   DummyInputStream(wigIn2)],
                                   mirrorScore = -1, mirror=True,
                                   debug = debug) :
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
