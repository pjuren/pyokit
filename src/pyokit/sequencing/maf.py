#!/usr/bin/python

"""
  Date of Creation: 23th May 2012
  Description:      Classes and functions for manipulating MAF files and
                    directories.

  Copyright (C) 2012-2014
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

import sys, os, random
from pyokit.io.bedIterators import BEDIterator
from pyokit.datastruct.intervalTree import IntervalTree
from pyokit.interface.cli import CLI, Option

DEFAULT_VERBOSITY = False

class MafLine :
  def __init__(self, line):
    parts = line.split()
    self.organism = parts[1].split(".")[0]
    self.chrom = ".".join(parts[1].split(".")[1:])
    self.start = int(parts[2])
    self.size = int(parts[3])
    self.strand = parts[4]
    self.val = parts[5]
    self.sequence = parts[6]

class MafBlock :
  """
    @summary: A <MafBlock> is a collection of <MafLine>s with a score. See
              constructor for concrete syntax example
  """

  def __init__(self, maflines, score, fh = None):
    """
      @summary: MafBlock constructor. Can build from a list of MafLines with
                a score, or can build by reading from a filehandle (or similar
                object). Concrete syntax example is:
                ... @todo: add syntax example
    """
    self.maflines = []
    self.score = None
    if fh != None :
      # ignore other things, load from stream
      line = None
      while True :
        line = fh.readline().strip()
        if line == "" : break
        self.maflines.append(MafLine(line))

      blockStarts_hg18 = [l.start for l in self.maflines if l.organism == "hg18"]
      blockEnds_hg18 = [l.start + l.size for l in self.maflines if l.organism == "hg18"]
      chroms = [l.chrom for l in self.maflines if l.organism == "hg18"]
      assert(len(blockStarts_hg18) == 1)
      assert(len(blockEnds_hg18) == 1)
      assert(len(chroms) == 1)
      self.start = blockStarts_hg18[0]
      self.end = blockEnds_hg18[0]
      self.chrom = chroms[0]
      #sys.stderr.write("loaded block for " + self.chrom + " " + str(self.start) + " " + str(self.end) + "\n")
    else :
      if score == None : raise ValueError("cannot make block with no score!")
      if len(maflines) <= 0 : raise ValueError("cannot make maf block with no lines!")
      self.maflines = maflines
      self.score = score

  def getLine(self, organism):
    """
      @summary: get the MafLine in this block which corresponds to the
                named organism
    """
    for line in self.maflines :
      if line.organism == organism : return line
    raise ValueError("no line for organism: " + organism)

  def getColumn(self, genomicCoord):
    """
      @return:  a list with the nucleotide for each organism at the given
                genomic coordinate.
    """
    human = [line for line in self.maflines if line.organism == "hg18"][0]

    relativeUngappedCoord = genomicCoord - human.start
    relativeGappedCoord = None
    nonGaps = 0
    for i in range(0, len(human.sequence)) :
      if human.sequence[i] != "-" : nonGaps += 1
      relativeGappedCoord = i
      if nonGaps == relativeUngappedCoord + 1 : break

    return [l.sequence[relativeGappedCoord] for l in self.maflines]

  def getColumnAsDictionary(self, genomicCoord):
    """
      @return: a dictionary indexed by organism name with the
               value of the nucleotide for each organism at the given
               genomic coordinate.
    """
    nucs = self.getColumn(genomicCoord)
    keys = [l.organism for l in self.maflines]
    return dict(zip(keys, nucs))


  def __str__(self):
    organismColWidth = max([len(x.organism + "." + x.chrom) for x in self.maflines]) + 4
    strandWidth = max([len(x.strand + " " + x.val) for x in self.maflines])
    startWidth = max([len(x.start) for x in self.maflines])
    sizeWidth = max([len(str(x.size)) for x in self.maflines])

    res = "a score=" + self.score + "\n"
    for line in self.maflines :
      paddingAmount = organismColWidth - len(str(line.organism + "." + line.chrom))
      padding = paddingAmount * " "
      strandPadAmount = strandWidth - len(line.strand + " " + line.val)
      strandPadding = strandPadAmount * " "
      orgToStartPadding = 4 * " "

      startPadding = " " * (startWidth - len(str(line.start)))
      sizePadding = " " * (sizeWidth - len(str(line.size)))

      #sys.stderr.write("before padding: " + line.organism + "." + line.chrom + str(line.start) + "\n")
      #sys.stderr.write("after  padding: " + line.organism + "." + line.chrom + padding + str(line.start) + "\n")
      #raw_input()
      res += "s " + line.organism + "." + line.chrom + padding + startPadding + str(line.start) + " " + sizePadding + str(line.size) + " " + str(line.strand) + strandPadding + " " + str(line.val) + " " + line.sequence + "\n"
    return res

class Maf :
  def __init__(self, fn):
    self.blocks = []
    self.loadMaf(fn)

  def __str__(self):
    return "\n".join([str(b) for b in self.blocks])

  def loadMaf(self, fn):
    readingBlock = False
    maflines = []

    for line in open(fn) :
      line = line.strip()

      # blank line marks end of block
      if line == "" :
        if readingBlock :
          self.blocks.append(MafBlock(maflines, blockScore))
          maflines = []
          blockScore = None
        readingBlock = False
      else :
        type = line.split()[0].strip()

        # an 'a' mark the start of a new block
        if type == "a" :
          if readingBlock == True : raise ValueError("found 'a' in the middle of a block")
          readingBlock = True
          blockScore = line.split("=")[1]

        # an 's' marks the content of a block
        if type == "s" : maflines.append(MafLine(line))

class MafIndex :
  def __init__(self, maf_fn, mafind_fn):
    self.indexElements = []
    self.maf_fn = maf_fn
    self.mafind_fn = mafind_fn
    self.loadIndex(self.mafind_fn)
    self.blockCache = None

  def loadIndex(self, fn):
    for e in BEDIterator(fn, scoreType = int) :
      self.indexElements.append(e)

  def getColumn(self, chrom, genomicCoord):
    # first see if it's a reference to what we already have loaded
    if self.blockCache == None or not (self.blockCache.chrom == chrom and genomicCoord >= self.blockCache.start and genomicCoord < self.blockCache.end) :
      indexLoc = None
      for e in self.indexElements :
        if genomicCoord >= e.start and genomicCoord <  e.end :
          indexLoc = e
          break
      if indexLoc == None : raise ValueError(str(genomicCoord) + " not in " + self.maf_fn)
      fh = open(self.maf_fn)
      fh.seek(indexLoc.score)
      self.blockCache = MafBlock(None, None, fh)
    return self.blockCache.getColumn(genomicCoord)

  def getColumnAsDictionary(self, chrom, genomicCoord):
    # first see if it's a reference to what we already have loaded
    if self.blockCache == None or not (self.blockCache.chrom == chrom and genomicCoord >= self.blockCache.start and genomicCoord < self.blockCache.end) :
      indexLoc = None
      for e in self.indexElements :
        if genomicCoord >= e.start and genomicCoord <  e.end :
          indexLoc = e
          break
      if indexLoc == None : raise ValueError(str(genomicCoord) + " not in " + self.maf_fn)
      fh = open(self.maf_fn)
      fh.seek(indexLoc.score)
      self.blockCache = MafBlock(None, None, fh)
    return self.blockCache.getColumnAsDictionary(genomicCoord)

class MafFile :
  def __init__(self, fn):
    self.fn = fn
    base, ext = os.path.splitext(os.path.split(fn)[-1])
    self.chrom = base.split(":")[0]
    self.start = int(base.split(":")[1].split("-")[0])
    self.end = int(base.split(":")[1].split("-")[1])

class MafDir :
  def __init__(self, dir, mafext, indext):
    maffiles = [os.path.join(dir,f) for f in os.listdir(dir) if os.path.splitext(f)[1] == "." + mafext]
    indexfiles = [os.path.join(dir,f) for f in os.listdir(dir) if os.path.splitext(f)[1] == "." + indext]

    self.mafFileTrees = {}
    mafFilesByChrom = {}
    for f in maffiles :
      mf = MafFile(f)
      if not mf.chrom in mafFilesByChrom : mafFilesByChrom[mf.chrom] = []
      mafFilesByChrom[mf.chrom].append(mf)
    for chrom in mafFilesByChrom :
      self.mafFileTrees[chrom] = IntervalTree(mafFilesByChrom[chrom])

    self.idxFileTrees = {}
    idxFilesByChrom = {}
    for f in indexfiles :
      idxf = MafFile(f)
      if not idxf.chrom in idxFilesByChrom : idxFilesByChrom[idxf.chrom] = []
      idxFilesByChrom[idxf.chrom].append(idxf)
    for chrom in idxFilesByChrom :
      self.idxFileTrees[chrom] = IntervalTree(idxFilesByChrom[chrom])

    self.mafIndexCached = None


  def getFilesFor(self, chrom, genomicCoord):
    if (not chrom in self.mafFileTrees) or (not chrom in self.idxFileTrees) : raise ValueError("Failed to find non-ambiguous files for " + chrom + " @ " + str(genomicCoord))
    mafFiles = self.mafFileTrees[chrom].intersectingPoint(genomicCoord)
    idxFiles = self.idxFileTrees[chrom].intersectingPoint(genomicCoord)
    if len(idxFiles) != 1 or len(mafFiles) != 1 : raise ValueError("Failed to find non-ambiguous files for " + chrom + " @ " + str(genomicCoord))
    return mafFiles[0].fn, idxFiles[0].fn

  def getColumn(self, chrom, genomicCoord):
    maf, ind = self.getFilesFor(chrom, genomicCoord)
    if self.mafIndexCached == None or not (self.mafIndexCached.maf_fn == maf and self.mafIndexCached.mafind_fn == ind) :
      #print "loading " + maf + " and " + ind
      self.mafIndexCached = MafIndex(maf, ind)
    return self.mafIndexCached.getColumn(chrom, genomicCoord)

  def getColumnAsDictionary(self, chrom, genomicCoord):
    maf, ind = self.getFilesFor(chrom, genomicCoord)
    if self.mafIndexCached == None or not (self.mafIndexCached.maf_fn == maf and self.mafIndexCached.mafind_fn == ind) :
      #print "loading " + maf + " and " + ind
      self.mafIndexCached = MafIndex(maf, ind)
    return self.mafIndexCached.getColumnAsDictionary(chrom, genomicCoord)
