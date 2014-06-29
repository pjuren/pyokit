#!/usr/bin/python

"""
  Date of Creation: 17th November 2010
  Description:      This module defines a class for buffered reading of
                    Wig data. At present, this requires the data to be
                    provided in a file, since we might make multiple
                    passes of it.

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

import sys, os, unittest, math

from pyokit.mapping.wigIterators import wigIterator, ITERATOR_SORTED_START
from pyokit.util.fileUtils import linesInFile, openFD
from pyokit.testing.dummyfiles import DummyInputStream, DummyOutputStream
from pyokit.datastruct.genomicMap import GenomicMap, GenomicMapError
from pyokit.wig import WigElement

DEFAULT_BUFFER_SIZE = 1000000
DEFAULT_CHUNK_PORTION = 1.5
DEFAULT_CHUNK_SIZE = math.ceil(DEFAULT_BUFFER_SIZE / DEFAULT_CHUNK_PORTION)
DEFAULT_MISSING_VAL = 0

class DummyIterator :
  def __init__(self, it):
    self.iterator = it
    self.lastChrom = None
    self.lastPos = None

  def ahead(self, chrom, pos):
    """
      returns true if the item would lie further ahead
    """
    if self.lastChrom == None or self.lastPos == None : return False
    if self.lastChrom < chrom : return True
    if self.lastChrom > chrom : return False
    if self.lastPos < pos : return True
    return False

  def next(self):
    n = self.iterator.next()
    self.lastChrom = n.chrom
    self.lastPos = n.start
    return n

class CircWig :
  def __init__(self, fn, verbose = False):
    self.fn = fn
    self.verbose = verbose
    self.resetIterator()

  def resetIterator(self):
    it = wigIterator(self.fn, sortedby = ITERATOR_SORTED_START, verbose = self.verbose)
    self.dummy = DummyIterator(it)

  def setIterator(self, chrom, pos):
    """
      set iterator relative to desired pos -- if it's behind, we
      make a new iterator at the start of the file, if it's ahead we
      return the iterator that is still unfinished
    """
    print "setting iterator relative to " + str(chrom) + "," + str(pos)
    if not self.dummy.ahead(chrom, pos) :
      print "desired spot is not ahead of last spot: " + str(self.dummy.lastChrom) + ", " + str(self.dummy.lastPos)
      self.resetIterator()

  def __iter__(self) :
    return self.dummy




class BufferedWig:
  def __init__(self, fn, bufferSize = DEFAULT_BUFFER_SIZE,
               chunkSize = DEFAULT_CHUNK_SIZE, verbose = False, debug = False):
    """
      @summary: constructor for a BufferedWig object
      @param fn: file descriptor for the underlying wig data.
                 can't be a stream because multiple passes might be made of
                 this data. Filename is fine.
      @param bufferSize: the max number of elements to keep in memory at any
                         given time
      @param chunkSize: the number of elements to flush from the buffer when
                        more space is needed.
      @param debug: output debugging messages if set to True
    """
    self.debug = debug
    if bufferSize != DEFAULT_BUFFER_SIZE and chunkSize == DEFAULT_CHUNK_SIZE:
      chunkSize = math.ceil(bufferSize / float(DEFAULT_CHUNK_PORTION))

    # debugging info
    self.misses = 0

    self.verbose = verbose

    self.missingVal = DEFAULT_MISSING_VAL
    self.chunkSize = chunkSize
    self.bufferSize = bufferSize
    self.fn = fn
    self.vals = {}
    self.buffer = GenomicMap()

    self.bufferMin = (None, None)
    self.bufferMax = (None, None)

    self.circwig = CircWig(self.fn, self.verbose)

    print "buffer is " + str(self.bufferSize)
    print "chunk size: " + str(self.chunkSize)

  def getItem(self, chrom, loc):
    """
      @summary: get the value of the item on chromosome <chrom> at location
                <loc> -- if it's not in the buffer, it'll be loaded first,
                if it's not in the wig, <self.missingVal> will be returned
      @note: the item is assumed to be only 1 position in length, with and
             end at start + 1
      @param chrom: the chromosome that the desired item is on
      @param loc: the genomic coordinate where the desired item starts
    """
    try :
      return self.buffer.getValue(chrom, loc)
    except GenomicMapError :
      if self.debug :
        des = str(chrom) + ", " + str(loc)
        longdes = des + " in \n" #+ str(self.buffer)
        #sys.stderr.write("failed to find " + longdes + "\n")

      # need to figure out whether we missed on this item because it
      # just isn't in the buffer (but might be in the wig elsewhere) or
      # because it's not in the wig at all...
      if self._inBuffer(chrom, loc) :
        if self.debug :
          pass
          #sys.stderr.write(des + " should be in buffer ") #+ str(self.buffer) + "\n")
        return WigElement(chrom, loc, loc+1, self.missingVal)
      self._loadItem(chrom, loc)
      try :
        if self.debug :
          longdes = des + " in \n" #+ str(self.buffer)
          sys.stderr.write("trying again for " + longdes + "\n")
        r = self.buffer.getValue(chrom, loc)
        if self.debug :
          longdes = des + " in \n" #+ str(self.buffer)
          sys.stderr.write("success for " + longdes + "\n")
        return r
      except GenomicMapError :
        if self.debug : sys.stderr.write("failed again, there is no entry in the wig for " + des + "\n")
        return WigElement(chrom, loc, loc+1, self.missingVal)

  def _inBuffer(self, chrom, loc):
    minChrom, minPos = self.bufferMin
    maxChrom, maxPos = self.bufferMax
    if minChrom == None or maxChrom == None or minPos == None or maxPos == None:
      return False
    if chrom < minChrom or chrom > maxChrom : return False
    if chrom == minChrom and loc < minPos : return False
    if chrom == maxChrom and loc > maxPos : return False
    return True

  def _loadItem(self, chrom, loc):
    """
      @summary: load the item on <chrom> at location <loc> into the buffer;
                Other surrounding items will also be loaded until the buffer
                has been filled, or the end of the stream is met.
      @param chrom: the chromosome that the desired item is on
      @param loc: the genomic coordinate where the desired item starts
      @note: the item is assumed to be just 1 position long, with an
             end at start + 1
    """
    self.misses += 1
    self.circwig.setIterator(chrom, loc)
    self.bufferMin = (chrom, loc)

    # if the buffer is full, make some space
    if len(self.buffer) >= self.bufferSize :
      if self.debug :
        sys.stderr.write("flushing buffer.. removing " +\
                         str(self.chunkSize) + " items \n")
      self.buffer.flush(self.chunkSize)
      if self.debug :
        sys.stderr.write("buffer now has " + str(len(self.buffer)) + " items")

    # find the item we need then keep loading items until we fill up the
    # buffer again or we run out of items
    if self.debug :
      sys.stderr.write("loading new data... ")
    adding = False
    first = True
    for item in self.circwig :
      if len(self.buffer) >= self.bufferSize : break
      if adding == False and item.chrom == chrom and item.start >= loc :
        adding = True
      if adding and not self.buffer.hasValue(item.chrom, item.start):
        self.buffer.addValue(item, item.chrom, item.start)
      first = False
    self.bufferMax = (item.chrom, item.start)
    if self.debug :
      sys.stderr.write("buffer now has " + str(len(self.buffer)) + " items")


class BufferedWigUnitTests(unittest.TestCase):
  """
    Unit tests for deadzone adjust
  """

  def setUp(self):
    # 59424 and 59427 are missing
    self.input = "chr1" + "\t" + "59420" + "\t" + "59421" + "\t" + "1" + "\n" +\
                 "chr1" + "\t" + "59421" + "\t" + "59422" + "\t" + "2" + "\n" +\
                 "chr1" + "\t" + "59422" + "\t" + "59423" + "\t" + "3" + "\n" +\
                 "chr1" + "\t" + "59423" + "\t" + "59424" + "\t" + "4" + "\n" +\
                 "chr1" + "\t" + "59425" + "\t" + "59426" + "\t" + "6" + "\n" +\
                 "chr1" + "\t" + "59426" + "\t" + "59427" + "\t" + "7" + "\n" +\
                 "chr1" + "\t" + "59428" + "\t" + "59429" + "\t" + "9" + "\n" +\
                 "chr1" + "\t" + "59429" + "\t" + "59430" + "\t" + "0" + "\n" +\
                 "chr2" + "\t" + "10000" + "\t" + "10001" + "\t" + "0" + "\n"

  def testInBuffer(self):
    """
      the desired item is in the buffer
    """
    infh = DummyInputStream(self.input)
    bw = BufferedWig(infh, bufferSize = 2, debug = False)
    self.assertTrue(bw.getItem("chr1", 59420).score == 1)
    self.assertTrue(bw.getItem("chr1", 59421).score == 2)
    self.assertTrue(bw.misses == 1)

  def testInWigNotInBuffer(self):
    """
      the desired item is in the wig, but hasn't been loaded yet
    """
    infh = DummyInputStream(self.input)
    bw = BufferedWig(infh, bufferSize = 2, debug = False)
    self.assertTrue(bw.getItem("chr1", 59422).score == 3)
    self.assertTrue(bw.misses == 1)
    self.assertTrue(bw.getItem("chr1", 59429).score == 0)
    self.assertTrue(bw.misses == 2)

  def testNotInWig(self):
    """
      the desired item does not exist in the wig at all
    """
    infh = DummyInputStream(self.input)
    bw = BufferedWig(infh, bufferSize = 2, debug = False)

    # case 1 -- it would have been in the buffer, so we shouldn't go looking for it
    bw.getItem("chr1", 59423)
    self.assertTrue(bw.getItem("chr1", 59424).score == DEFAULT_MISSING_VAL)
    self.assertTrue(bw.misses == 1)

    # case 2 -- it wouldn't be in the buffer, so we go looking
    self.assertTrue(bw.getItem("chr1", 59427).score == DEFAULT_MISSING_VAL)
    self.assertTrue(bw.misses == 2)

  def testAfterEndOfChrom(self):
    """
      asking for an item beyond the end of the chrom should not cause
      reloads
    """
    infh = DummyInputStream(self.input)
    bw = BufferedWig(infh, bufferSize = 2, debug = False)
    bw.getItem("chr1", 59429)
    self.assertTrue(bw.getItem("chr1", 59490).score == DEFAULT_MISSING_VAL)
    self.assertTrue(bw.getItem("chr2", 10000).score == 0)
    self.assertTrue(bw.misses == 1)

  def testBeforeFirstItem(self):
    """
      asking for an item that falls before the first item in the wig
      should cause miss the first time, but not again
    """
    infh = DummyInputStream(self.input)
    bw = BufferedWig(infh, bufferSize = 2, debug = False)
    self.assertTrue(bw.getItem("chr1", 10).score == DEFAULT_MISSING_VAL)
    self.assertTrue(bw.misses == 1)
    self.assertTrue(bw.getItem("chr1", 90).score == DEFAULT_MISSING_VAL)
    self.assertTrue(bw.misses == 1)

  def testPassedVal(self):
    """
      asking for a value that we already passed
    """
    infh = DummyInputStream(self.input)
    bw = BufferedWig(infh, bufferSize = 2, debug = False)
    bw.getItem("chr1", 59428)
    self.assertTrue(bw.getItem("chr1", 59426).score == 7)
    self.assertTrue(bw.getItem("chr2", 10000).score == 0)
    self.assertTrue(bw.misses == 3)

if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])
