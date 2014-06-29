#!/usr/bin/python
"""
  Date of Creation: 15th Oct 2011
  Description:      Class for loading a whole Wig file at once and
                    performing random access to the elements.
                    Lookup in O(log(n)) for a single element
                    Building data structure is O(nlog(n)).

  Copyright (C) 2011-2014
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

import sys, unittest
from pyokit.io.wigIterators import wigIterator
from pyokit.datastruct.intervalTree import IntervalTree
from pyokit.testing.dummyfiles import DummyInputStream, DummyOutputStream

class WigFileError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)

class WigFile :
  def __init__(self, filename, verbose = False):
    """
      @summary: constructor for the WigFile class
      @param filename: can be either a string filename or a
                       file-handle-like object
    """
    self.filename = filename
    self.itrees = {}
    self.verbose = verbose
    self.__load(verbose=self.verbose)

  def __load(self, verbose=False):
    """
      @summary: load the contents of a wig file into this object
    """
    byChrom = {}
    for e in wigIterator(self.filename, verbose=verbose) :
      if not e.chrom in byChrom : byChrom[e.chrom] = []
      byChrom[e.chrom].append(e)
    for chrom in byChrom :
      self.itrees[chrom] = IntervalTree(byChrom[chrom], openEnded=True)

  def contains(self, chrom, point):
    """
      @summary: return True if this WigFile has an element covering the
                given position, False otherwise
    """
    if not chrom in self.itrees : return False
    hits = self.itrees[chrom].intersectingPoint(point)
    return len(hits) > 0

  def getElement(self, chrom, point):
    """
      @summary: get the WigElement that is intersected by the given point
                returns None if there is no such element
      @return:  a WigElement
      @raise WigFileError: if more than one element intersects the point
    """
    if not chrom in self.itrees : return None
    hits = self.itrees[chrom].intersectingPoint(point)
    if len(hits) == 0 : return None
    if len(hits) > 1  : raise WigFileError("multiple entries intersect " +\
                                           str(chrom) + " at " + str(point))
    return hits[0]

  def getScore(self, chrom, point):
    """
      @summary: get the value (float) that is intersected by the given point
                returns None if there is no such element
      @return:  a float
      @raise WigFileError: if more than one element intersects the point
    """
    e = self.getElement(chrom, point)
    if e == None : return None
    return e.score


class WigFileUnitTests(unittest.TestCase):
  """
    @summary: Unit tests for WigFile
  """

  def setUp(self):
    pass

  def testGetValue(self):
    debug = False
    wigIn = "chr1" + "\t" + "01" "\t" + "10" + "\t" + "5" + "\n" +\
            "chr1" + "\t" + "40" "\t" + "50" + "\t" + "3" + "\n" +\
            "chr2" + "\t" + "30" "\t" + "34" + "\t" + "2" + "\n" +\
            "chr4" + "\t" + "30" "\t" + "60" + "\t" + "6" + "\n"
    test = [("chr1", 4,5),
            ("chr2",31,2),
            ("chr4",40,6),
            ("chr2",30,2),
            ("chr1",45,3)]

    infh = DummyInputStream(wigIn)
    wf = WigFile(infh)
    for chrom, point, ans in test :
      self.assertTrue(wf.getScore(chrom, point) == ans)

  def testFailOverlap(self):
    debug = False
    wigIn = "chr1" + "\t" + "01" "\t" + "10" + "\t" + "5" + "\n" +\
            "chr1" + "\t" + "09" "\t" + "50" + "\t" + "3" + "\n"
    infh = DummyInputStream(wigIn)
    wf = WigFile(infh)
    self.assertRaises(WigFileError, wf.getScore, "chr1", 9)

if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])
