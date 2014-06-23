#!/usr/bin/python

"""
  Date of Creation: 24th November 2010
  Description:      Classes and functions for manipulating refseq
                    transcriptome references

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

import sys, os, unittest
from pyokit.datastruct.intervalTree import IntervalTree
from pyokit.mapping.bedIterators import BEDIterator
from pyokit.mapping.bed import toGenomicCoordinates, BEDElement
from pyokit.mapping.transcript import Transcript, TranscriptError
from pyokit.testing.dummyfiles import DummyInputStream, DummyOutputStream

KEY_SEP = "___"

class ReferenceMap :
  """
    @summary: A ReferenceMap is a collection of transcripts arranged as a Tree
              so that transcripts overlapping regions can be found quickly
  """

  def __init__(self, transcripts):
    """
      @summary: construct a transcript map from the provided list of
                transcripts
    """
    self.trees = {}
    chroms = set([t.chrom for t in transcripts])
    for chrom in chroms :
      ts = [t for t in transcripts if t.chrom == chrom]
      self.trees[chrom] = IntervalTree(ts)


  def getTranscriptsAt(self, chrom, start, end = None):
    """
      @summary: get all transcripts that overlap the given point or interval
      @return: A list of transcripts overlapping the region
    """
    tree = self.trees[chrom]
    if end == None : return tree.intersectingPoint(start)
    else : return tree.intersectingInterval(start,end)

  def getRelativeIndices(self, chrom, pos, debug = False):
    """
      @summary: given a genomic coordinate, determine the relative coordinate
                from the start of the transcript it's in. If it overlaps
                multiple transcripts, a list of relative coordinates will be
                returned. If the coordinate overlaps no transcripts, None will
                be returned
      @return: a list of relative positions within each transcript that are
               overlapped by this position
      @raise ReferenceError: if the position given is not in any exon
                             contained in the ReferenceMap
    """
    res = []
    transcripts = self.getTranscriptsAt(chrom, pos)
    if transcripts == None :
      raise ReferenceError(chrom + " at " + str(pos) + " " +\
                           "is not in the transcriptome represented " +\
                           "by this ReferenceMap -- no transcripts")
    found = False
    for t in transcripts :
      try :
        res.append(t.absoluteToRelative(chrom, pos, debug = debug))
        found = True
      except TranscriptError :
        pass

    if not found :
      sys.stderr.write("\nfailed to find matching exon in transcript hits as follows\n")
      for t in transcripts :
        sys.stderr.write(str(t) + "\n\n")
      raise ReferenceError(chrom + " at " + str(pos) + " " +\
                           "is not in the transcriptome represented " +\
                           "by this ReferenceMap -- no exons")
    return res



def readReferenceAsTranscripts(reffn, verbose = False, overlapWarn = False, strandWarn = False):
  """
    @param overlapWarn: if True, transcripts with overlapping exons are not
                        considered a fatal error and a warning is issued
                        instead of raising an exception
    @param strandWarn: if True, treat mixed strand in exon as non-fatal
                         and just issue a warning
  """
  transcripts = []
  ts = readTranscriptomeByGene(reffn, verbose = verbose)
  for k in ts :
    transcripts.append(Transcript(ts[k], overlapWarn, strandWarn))
  return transcripts


def makeReferenceKey(element):
  return element.chrom + KEY_SEP + element.name.split("_exon_")[0]

def getUniqueTranscriptIDs(elements):
  """
    @summary: given a list of exons taken from a transcriptome reference,
              get the list of chrom and refseq keys that uniquely ID
              each element
  """
  unique = set([element.chrom + KEY_SEP + element.name.split("_exon_")[0]
                for element in elements])
  return list(unique)

def readTranscriptomeByGene(reffn, verbose):
  """
    @summary: read a transcriptome reference and group elements by gene
    @return: returns a dictionary where each item is a list of exons
             for a given gene. Entries are indexed by keys which are made
             by concatenating the chromosome, refseq and strand separated by
             KEY_SEP.
  """
  t = readTranscriptomeFlat(reffn, verbose = verbose)
  return splitReferenceByGene(t)

def splitReferenceByGene(elements, verbose = False):
  """
    @summary: take a list of elements from a transcritome reference
              and split into lists for each refseq name, collected
              into a dictionary. Entries are indexed by keys which are made
              by concatenating the chromosome, refseq and strand separated by
              KEY_SEP.
  """
  transcripts = {}
  for element in elements :
    refseq = element.name.split("_exon_")[0]
    chrom = element.chrom
    strand = element.strand.strip()
    key = chrom + KEY_SEP + refseq + KEY_SEP + strand

    if not key in transcripts : transcripts[key] = []
    transcripts[key].append(element)
  return transcripts

def readTranscriptomeFlat(reffn, verbose = False):
  """
    @summary: read a transcriptome reference and return all elements
              as a flat list
    @return: list of all elements
  """
  return [element for element in BEDIterator(reffn, verbose = verbose)]




class ReferenceUnitTests(unittest.TestCase):
  """
    Unit tests for Reference
  """

  def setUp(self):
    self.rawinput = "chr1" + "\t" + "10" + "\t" + "20" + "\t" + "gene1_exon_1" + "\t" + "0" + "\t" + "+" + "\n" +\
                    "chr1" + "\t" + "40" + "\t" + "50" + "\t" + "gene1_exon_2" + "\t" + "0" + "\t" + "+" + "\n" +\
                    "chr1" + "\t" + "60" + "\t" + "70" + "\t" + "gene1_exon_3" + "\t" + "0" + "\t" + "+" + "\n" +\
                    "chr2" + "\t" + "10" + "\t" + "20" + "\t" + "gene2_exon_1" + "\t" + "0" + "\t" + "+" + "\n" +\
                    "chr2" + "\t" + "30" + "\t" + "40" + "\t" + "gene2_exon_2" + "\t" + "0" + "\t" + "+" + "\n" +\
                    "chr2" + "\t" + "50" + "\t" + "60" + "\t" + "gene2_exon_3" + "\t" + "0" + "\t" + "+" + "\n" +\
                    "chr2" + "\t" + "70" + "\t" + "80" + "\t" + "gene2_exon_4" + "\t" + "0" + "\t" + "+" + "\n" +\
                    "chr2" + "\t" + "110" + "\t" + "120" + "\t" + "gene2_exon_1" + "\t" + "0" + "\t" + "-" + "\n" +\
                    "chr2" + "\t" + "130" + "\t" + "140" + "\t" + "gene2_exon_2" + "\t" + "0" + "\t" + "-" + "\n" +\
                    "chr2" + "\t" + "150" + "\t" + "160" + "\t" + "gene2_exon_3" + "\t" + "0" + "\t" + "-" + "\n" +\
                    "chr2" + "\t" + "170" + "\t" + "180" + "\t" + "gene2_exon_4" + "\t" + "0" + "\t" + "-" + "\n" +\
                    "chr2" + "\t" + "45" + "\t" + "55" + "\t" + "gene3_exon_1" + "\t" + "0" + "\t" + "+" + "\n" +\
                    "chr2" + "\t" + "90" + "\t" + "95" + "\t" + "gene3_exon_2" + "\t" + "0" + "\t" + "+" + "\n"
    #-----

    e1 = BEDElement("chr1", 10, 20, "gene1_exon_1", 0, "+")
    e2 = BEDElement("chr1", 40, 50, "gene1_exon_2", 0, "+")
    e3 = BEDElement("chr1", 60, 70, "gene1_exon_3", 0, "+")

    e4 = BEDElement("chr2", 10, 20, "gene2_exon_1", 0, "+")
    e5 = BEDElement("chr2", 30, 40, "gene2_exon_2", 0, "+")
    e6 = BEDElement("chr2", 50, 60, "gene2_exon_3", 0, "+")
    e7 = BEDElement("chr2", 70, 80, "gene2_exon_4", 0, "+")
    e41 = BEDElement("chr2", 110, 120, "gene2_exon_1", 0, "-")
    e51 = BEDElement("chr2", 130, 140, "gene2_exon_2", 0, "-")
    e61 = BEDElement("chr2", 150, 160, "gene2_exon_3", 0, "-")
    e71 = BEDElement("chr2", 170, 180, "gene2_exon_4", 0, "-")

    e8 = BEDElement("chr2", 45, 55, "gene3_exon_1", 0, "+")
    e9 = BEDElement("chr2", 90, 95, "gene3_exon_2", 0, "+")

    self.t1 = Transcript([e1, e2, e3])
    self.t2 = Transcript([e4, e5, e6, e7])
    self.t3 = Transcript([e8, e9])
    self.t4 = Transcript([e41,e51,e61,e71])
    self.r1 = ReferenceMap([self.t1, self.t2, self.t3])

  def testReadReferenceAsTranscripts(self):
    debug = True
    transcripts = readReferenceAsTranscripts(DummyInputStream(self.rawinput))
    expect = [self.t1,self.t2,self.t3,self.t4]
    transcripts.sort(key = lambda x: x.transcriptID())
    expect.sort(key = lambda x: x.transcriptID())
    if debug :
      sys.stderr.write("\nwe expected to get.. \n ----------- \n")
      for t in expect :
        sys.stderr.write(str(t) + "\n")
      sys.stderr.write("\nwe got.. \n ----------- \n")
      for t in transcripts :
        sys.stderr.write(str(t) + "\n")
    self.assertTrue(expect == transcripts)


  def testGetRelativeIndices(self):
    debug = False
    fail = None

    for chrom, point, res in [("chr1", 15, [5]),
                              ("chr1", 25, fail),
                              ("chr2", 51, [21,6])] :
      if res == fail :
        self.assertRaises(ReferenceError, self.r1.getRelativeIndices,
                          chrom, point, debug)
      else :
        rel = self.r1.getRelativeIndices(chrom, point, debug = debug)
        if debug :
          sys.stderr.write("rel for abs " + chrom + " at " + str(point) +\
                           " is " + str(rel) + " and we're expecting " +\
                           str(res) + "\n")
        self.assertTrue(rel == res)


if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])
