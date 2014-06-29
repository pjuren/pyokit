#!/usr/bin/python
"""
  Date of Creation: 1st Oct 2011
  Description:      An IndexedWig allows fast random access to a large Wig
                    file without requiring the whole file to be loaded into
                    memory at any given time. It builds an index of the wig
                    file and uses a cache to store recently accessed blocks.
                    The file must be sorted. Performance is reasonable under
                    the assumption of locality of reference (next reference
                    is expected to be near current one), but complete random
                    access could force a load on each lookup, causing poor
                    performance.

  Copyright (C) 2011-2014
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

import sys, math, os, unittest, time
from random import shuffle
from heapq import *
from pyokit.testing.dummyfiles import DummyInputStream, DummyOutputStream
from pyokit.datastruct.genomicInterval import parseWigString
from pyokit.datastruct.intervalTree import IntervalTree
from pyokit.util.progressIndicator import ProgressIndicator

class IndexedWigError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)

class WigBlock :
  """
    everything in the block is on the same chrom
    it includes all elements between start (inclusive) and end (exclusive)
  """

  def __init__(self, fileLocation, e, blocksize, debug = False):
    self.fileLocation = fileLocation
    self.capacity = blocksize
    self.chrom = e.chrom
    self.start = e.start
    self.debug = debug
    self.iTree = None
    self.end = e.end
    self.data = []
    self.size = 0


  def __len__(self):
    """
      @summary: return the size of this block. The size is defined as the
                total number of wig elements that would be in this block
                regardless of whether they are loaded now or not
    """
    return self.size

  def isfull(self):
    return len(self) >= self.capacity

  def isempty(self):
    return len(self) == 0

  def add(self, e):
    if e.chrom != self.chrom :
      msg = "trying to add item on chrom " + str(e.chrom) + "to WigBlock " +\
            "for " + str(self.chrom)
      raise IndexedWigError(msg)
    self.end = e.end
    self.size += 1

  def lookup(self, chrom, point):
    if chrom != self.chrom :
      msg = "looking for value on chrom " + chrom + " in " + str(self) +\
            " failed; chromosomes don't match"
      raise IndexedWigError(msg)
    hits = self.iTree.intersectingPoint(point)
    if len(hits) == 0 :
      msg = "lookup for value on chrom " + str(chrom) + " at location " +\
             str(point) + " failed; location not covered in wig file"
      raise IndexedWigError(msg)
    elif len(hits) > 1 :
      msg = "lookup for value on chrom " + str(chrom) + " at location " +\
             str(point) + " failed; multiple entries intersecting this point"
      raise IndexedWigError(msg)
    return hits[0]

  def contains(self, chrom, point):
    if chrom != self.chrom : return False
    hits = self.iTree.intersectingPoint(point)
    if len(hits) == 0 : return False
    elif len(hits) > 1 :
      msg = "lookup for value on chrom " + str(chrom) + " at location " +\
             str(point) + " failed; multiple entries intersecting this point" +\
             "\n".join([str(e) for e in hits])
      raise IndexedWigError(msg)
    return True

  def isPopulated(self):
    return len(self.data) != 0

  def populate(self, filehandle):
    #print "populating " + str(self)
    #print "seeking to " + str(self.fileLocation)
    filehandle.seek(self.fileLocation)
    if self.debug : sys.stderr.write("populating " + str(self) + "\n")
    for line in filehandle :
      # get next element
      line = line.strip()
      if self.debug :
        sys.stderr.write("\t" + "current line is " + str(line) + "\n")
      line = line.strip()
      if line == "": continue
      e = parseWigString(line)

      # we're done if we've left this block's chrom, or if we've moved beyond
      # the end of this blocks boundary.
      if e.chrom != self.chrom or e.start > self.end : break
      self.data.append(e)
    if self.debug : sys.stderr.write("built tree for " + str(self) + "\n")
    if len(self.data) == 0 : print "empty! --> " + str(self)
    self.iTree = IntervalTree(self.data, openEnded=True)

  def depopulate(self):
    self.data = []
    self.iTree = None

  def __str__(self):
    res = "Block at file location " + str(self.fileLocation) +\
          " has capacity " + str(self.capacity) + " is on chrom " +\
          str(self.chrom) + " starting at " + str(self.start) +\
          " and ending at " + str(self.end)
    return res


class IndexedWig :
  def __init__(self, filename, blocksize, blockCapacity,
               debug=False, verbose=False):
    """
      @summary: Constructor for IndexedWig objects
      @param filename: handle, or filename for the wig file
      @param blocksize: number of wig elements to store in each block
      @param blockCapacity: number of blocks to keep in memory at a time
      @param debug: if true, output additional debugging info
      @param verbose: if true, output additional run info when building object
    """
    self.itrees = {}
    self.blocksize = blocksize
    self.blockCapacity = blockCapacity
    self.populatedBlocks = []
    self.blocksByChrom = {}
    self.debug = debug
    self.verbose = verbose
    self.handle = filename
    if type(self.handle).__name__ == "str" :
      self.handle = open(self.handle)
    self.build()

  def build(self):
    currentBlock = None
    at = self.handle.tell()
    seenChroms = set()
    lastIndexSeen = -1

    if self.verbose :
      try :
        pind = ProgressIndicator(totalToDo = os.path.getsize(self.handle.name),
                                       messagePrefix = "completed",
                                       messageSuffix = "of building index for " +\
                                                        self.handle.name)
      except :
        sys.stderr.write("IndexedWig -- warning: " +\
                         "unable to show progress for stream\n")
        self.verbose = False

    ### note, for loop seems to buffer the file and so tell() gives a
    ### location that is not where the current line was read from, so
    ### we stick to readline instead.
    rline = None
    while rline != "" :
      # get the next element
      rline = self.handle.readline()
      line = rline.strip()
      if line == "": continue
      e = parseWigString(line)

      # keep track of what chroms we've seen for checking order
      if not e.chrom in seenChroms :
        seenChroms.add(e.chrom)
        lastIndexSeen = -1

      # check chrom order is ok
      for seenChrom in seenChroms :
        if seenChrom > e.chrom :
          msg = "wig file is not sorted, entry for chrom " + str(seenChrom) +\
                " appears after entry for " + str(e.chrom)
          raise IndexedWigError(msg)
      # check position order is ok
      if e.start < lastIndexSeen :
        msg = "wig file is not sorted, entry for chrom " + str(e.chrom) +\
                " at " + str(e.start) + " appears after " + str(lastIndexSeen)
        raise IndexedWigError(msg)

      # update the last index we've seen
      lastIndexSeen = e.end

      # debugging message if the current block is full
      if self.debug == True :
        sys.stderr.write("processing " + str(e))
        if currentBlock != None :
          sys.stderr.write("; is current block full?" +\
                           str(currentBlock.isfull()) + "\n")
        else :
          sys.stderr.write("\n")

      # we might need to make a new block for this element
      if currentBlock == None or currentBlock.isfull() or \
         currentBlock.chrom != e.chrom :
        if self.debug : sys.stderr.write("making new block with " + str(e) + "\n")
        if currentBlock != None :
          if self.debug : sys.stderr.write("closed block: " + str(currentBlock) + "\n")
          if not currentBlock.chrom in self.blocksByChrom :
            self.blocksByChrom[currentBlock.chrom] = []
          self.blocksByChrom[currentBlock.chrom].append(currentBlock)
        currentBlock = WigBlock(at, e, self.blocksize)

      # add the element to the current block
      currentBlock.add(e)

      at = self.handle.tell()

      if self.verbose :
        pind.done = self.handle.tell()
        pind.showProgress()

    # don't forget to add the final block
    if currentBlock != None :
      if self.debug : sys.stderr.write("closed block: " + str(currentBlock) + "\n")
      if not currentBlock.chrom in self.blocksByChrom :
        self.blocksByChrom[currentBlock.chrom] = []
      self.blocksByChrom[currentBlock.chrom].append(currentBlock)

    # build the interval trees
    for chrom in self.blocksByChrom :
      self.itrees[chrom] = IntervalTree(self.blocksByChrom[chrom], openEnded=True)


  def contains(self, chrom, point):
    """
      @summary: returns true if the indexed wig contains a value for the
                given location, false otherwise
    """
    if not chrom in self.itrees : return False
    blocks = self.itrees[chrom].intersectingPoint(point)
    if len(blocks) == 0 : return False

    if len(blocks) > 1 :
      msg = "lookup for value on chrom " + str(chrom) + " at location " +\
             str(point) + " failed; multiple entries intersecting this point"
      raise IndexedWigError(msg)
    block = blocks[0]

    if not block.isPopulated() :
      block.populate(self.handle)
      if len(self.populatedBlocks) >= self.blockCapacity :
        if self.debug : sys.stderr.write("popping old block off heap\n")
        t, o = heappop(self.populatedBlocks)
        o.depopulate()
      heappush(self.populatedBlocks, (time.time, block))

    return block.contains(chrom, point)


  def lookup(self, chrom, point):
    """
      @summary: get the wig element which intersects the given point
    """
    if not chrom in self.itrees :
      msg = "lookup for value on chrom " + str(chrom) + " failed; no such " +\
            "chromosome in wig file"
      raise IndexedWigError(msg)

    # get the block that this point is in
    blocks = self.itrees[chrom].intersectingPoint(point)
    if self.debug :
      sys.stderr.write("lookup for " + str(chrom) + " at " +\
                        str(point) + " got these blocks\n" +\
                       "\n".join([str(x) for x in blocks]) + "\n")
    if len(blocks) == 0 :
      msg = "lookup for value on chrom " + str(chrom) + " at location " +\
             str(point) + " failed; location not covered in wig file"
      raise IndexedWigError(msg)
    elif len(blocks) > 1 :
      msg = "lookup for value on chrom " + str(chrom) + " at location " +\
             str(point) + " failed; multiple entries intersecting this point"
      raise IndexedWigError(msg)
    block = blocks[0]

    #print "looking up " + str(chrom) + " at " + str(point) + " got us block " + str(block)

    # populate the block if it's not already; age out old blocks from
    # our cache if we need to make new space
    if self.debug :
      sys.stderr.write("block " + str(block) + " is populated? " +\
                        str(block.isPopulated()) + "\n")
    if not block.isPopulated() :
      block.populate(self.handle)
      if len(self.populatedBlocks) >= self.blockCapacity :
        if self.debug : sys.stderr.write("popping old block off heap\n")
        t, o = heappop(self.populatedBlocks)
        o.depopulate()
      heappush(self.populatedBlocks, (time.time, block))

    return block.lookup(chrom, point)

  def __str__(self):
    """
      @summary: produce a string representation of this IndexedWig object.
                We just show the details of the blocks, we don't show
                individual wig elements since these are often not in memory
    """
    res = ""
    for chrom in self.blocksByChrom :
      for block in self.blocksByChrom[chrom] :
        res += (str(block) + "\n")
    return res


class TestIndexedWig(unittest.TestCase):
  def setUp(self):
    self.input = ["chr1" + "\t" + "10" + "\t" + "11" + "\t" + "1"  + "\n",
                  "chr1" + "\t" + "30" + "\t" + "31" + "\t" + "2"  + "\n",
                  "chr1" + "\t" + "45" + "\t" + "46" + "\t" + "3"  + "\n",
                  "chr1" + "\t" + "96" + "\t" + "97" + "\t" + "4"  + "\n",
                  "chr2" + "\t" + "34" + "\t" + "35" + "\t" + "5"  + "\n",
                  "chr2" + "\t" + "46" + "\t" + "47" + "\t" + "6"  + "\n",
                  "chr2" + "\t" + "57" + "\t" + "58" + "\t" + "7"  + "\n",
                  "chr2" + "\t" + "69" + "\t" + "70" + "\t" + "8"  + "\n",
                  "chr3" + "\t" + "79" + "\t" + "80" + "\t" + "9"  + "\n",
                  "chr3" + "\t" + "83" + "\t" + "84" + "\t" + "10" + "\n",
                  "chr3" + "\t" + "92" + "\t" + "93" + "\t" + "11" + "\n"]

  def testRecall(self) :
    """
      @summary: here we're just testing that we can get back what we put in
                to the IndexedWig object in a random order
    """
    debug = False
    infh = DummyInputStream("\n".join(self.input))
    answers = [parseWigString(l) for l in self.input if l.strip() != ""]
    shuffle(answers)

    iwig = IndexedWig(infh, 2, 2, debug, verbose=False)
    if debug : sys.stderr.write("iwig structure is: \n" + str(iwig) + "\n")
    print "done"
    for e in answers :
      ans = iwig.lookup(e.chrom, e.start)
      if debug :
        sys.stderr.write("expect: " + str(e.score) + ", got: " + str(ans) + "\n")
      assert(e.score == ans.score)

  def testOutOfOrder(self):
    """
      @summary: make sure the right exception is thrown if the wig file is not
                sorted
    """
    shuffle(self.input)
    debug = False
    infh = DummyInputStream("\n".join(self.input))
    self.assertRaises(IndexedWigError, IndexedWig, infh, 2, 2, debug)

if __name__ == '__main__':
  unittest.main()
