#!/usr/bin/python
"""
  Date of Creation: 15th Oct 2011
  Description:      A WigDir is a directory containing a collection of
                    wig files split by chromosome

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

import sys, unittest, os
from pyokit.testing.dummyfiles import DummyInputStream, DummyOutputStream
from pyokit.io.wigFile import WigFile

class WigDir :
  def __init__(self, dir, extension="wig", verbose=False):
    self.extension = extension
    self.files = {}
    self.currentlyLoaded = None
    self.verbose=verbose
    self.__load(dir, verbose)

  def contains(self, chrom, point):
    """
      @summary: return True if this WigDir has an element covering the
                given position, False otherwise
    """
    filename = self._chromToFilename(chrom)
    if filename == None : return False
    if self.currentlyLoaded != filename : self._loadFile(filename)
    return self.currentlyLoaded.contains(chrom, point)

  def getElement(self, chrom, point):
    """
      @summary: get the WigElement that is intersected by the given point
                returns None if there is no such element
      @return:  a WigElement
      @raise WigFileError: if more than one element intersects the point
    """
    fh = self.__chromToFileHandle(chrom)
    if fh == None : return None
    if self.currentlyLoaded != fh : self.__loadFile(fh, verbose=self.verbose)
    return self.currentlyLoaded.getElement(chrom, point)

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

  def __str__(self):
    return ",".join(self.files.keys())

  def __chromToFileHandle(self, chrom):
    """
      @summary: find the file that contains entries for the given chromosome
    """
    if not chrom in self.files : return None
    return self.files[chrom]
    #hits = [f for f in self.files.keys() if f.name.find(chrom) != -1]
    #mhits = [f for f in hits if len(f.nam) == max(hits, key=len)]
    #if len(mhits) == 0 : return None
    #if len(mhits) > 1  : return None
    #return mhits[0]

  def __load(self, dir, verbose=False):
    """
      @summary: builds self.files, which is a dictionary of file-handle-like
                objects indexed by name.
      @param dir: -if single string we treat as directory of files to load
                  -if a list of string, we treat as list of filenames to load
                  -if list of anything else, treat as list of file-handle-like
                   objects (must have .name attributes)
    """
    if type(dir).__name__ == "str" :
      for f in os.listdir(dir) :
        name = os.path.splitext(f)[0]
        self.files[name] = open(os.path.join(dir,f))
    else :
      for f in dir :
        if type(f).__name__ == "str" : f = open(f)
        name = os.path.splitext(f.name)[0]
        self.files[name] = f

  def __loadFile(self, filename, verbose=False):
    """
      @summary: load the specified file.
      @param filename: if a string, we treat it as a name and try to open a
                       handle to it. Any other object is treated as file-like
                       and we try to seek the start of the file and build a
                       WigFile object from it
    """
    fh = None
    if type(filename).__name__ == "str" : fh = open(filename)
    else :
      fh = filename
      fh.seek(0)
    self.currentlyLoaded = WigFile(fh, verbose=verbose)


class WigDirUnitTests(unittest.TestCase):
  """
    @summary: Unit tests for WigFile
  """

  def setUp(self):
    pass

  def testGetScore(self):
    debug = False
    wig1 = DummyInputStream(
           "chr1" + "\t" + "01" "\t" + "10" + "\t" + "01" + "\n" +\
           "chr1" + "\t" + "40" "\t" + "50" + "\t" + "02" + "\n" +\
           "chr1" + "\t" + "55" "\t" + "62" + "\t" + "03" + "\n" +\
           "chr1" + "\t" + "74" "\t" + "87" + "\t" + "04" + "\n")
    wig2 = DummyInputStream(
           "chr2" + "\t" + "06" "\t" + "11" + "\t" + "05" + "\n" +\
           "chr2" + "\t" + "33" "\t" + "43" + "\t" + "06" + "\n" +\
           "chr2" + "\t" + "56" "\t" + "63" + "\t" + "07" + "\n" +\
           "chr2" + "\t" + "77" "\t" + "78" + "\t" + "08" + "\n")
    wig3 = DummyInputStream(
           "chr3" + "\t" + "03" "\t" + "09" + "\t" + "09" + "\n" +\
           "chr3" + "\t" + "15" "\t" + "16" + "\t" + "10" + "\n" +\
           "chr3" + "\t" + "22" "\t" + "25" + "\t" + "11" + "\n" +\
           "chr3" + "\t" + "38" "\t" + "39" + "\t" + "12" + "\n")
    wig4 = DummyInputStream(
           "chr4" + "\t" + "11" "\t" + "12" + "\t" + "13" + "\n" +\
           "chr4" + "\t" + "12" "\t" + "13" + "\t" + "14" + "\n" +\
           "chr4" + "\t" + "13" "\t" + "14" + "\t" + "15" + "\n" +\
           "chr4" + "\t" + "14" "\t" + "15" + "\t" + "16" + "\n")
    wig1.name="chr1.wig"
    wig2.name="chr2.wig"
    wig3.name="chr3.wig"
    wig4.name="chr4.wig"
    wigInfs = [wig1, wig2, wig3, wig4]
    test = [("chr1", 4, 1),
            ("chr1",45, 2),
            ("chr2",34, 6),
            ("chr2",77, 8),
            ("chr4",14,16)]

    wd = WigDir(wigInfs)
    for chrom, point, ans in test :
      self.assertTrue(wd.getScore(chrom, point) == ans)

if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])
