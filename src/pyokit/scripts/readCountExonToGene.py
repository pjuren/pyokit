#!/usr/bin/python

"""
  Date of Creation: 25rd Jun 2011

  Description:    Given a tab delineated file with two coloumns, one for exon
                  name and the other for read count, and a BED file indicating
                  the locations of genes, produce gene read counts

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
"""

# standard python imports
import sys
import os
import unittest

# pyokit imports
from pyokit.interface.cli import CLI, Option
from pyokit.testing.dummyfiles import DummyInputStream, DummyOutputStream
from pyokit.io.bedIterators import BEDIterator

###############################################################################
#                                  CONSTANTS                                  #
###############################################################################

DEFAULT_VERBOSITY = False


###############################################################################
#                              EXCEPTION CLASSES                              #
###############################################################################

class readCountExonToGeneError(Exception):
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


###############################################################################
#                                 EXON CLASS                                  #
###############################################################################

class Exon:
  def __init__(self, details, count):
    # NM_000014_exon_16_0_chr12_9134219_r
    self.details = details
    head, tail = details.split("_exon_")
    self.refseq = head.strip()
    if tail.split("_")[-1] == "r":
      self.strand = "-"
    else:
      self.strand = "+"
    self.exonNumber = int(tail.split("_")[0])
    self.loc = int(tail.split("_")[-2])
    self.chrom = "_".join(tail.split("_")[2:-2]).strip()
    self.count = int(count)

    # sanity check
    if not (str(self) == str(self.details + "\t" + str(count))):
      raise readCountExonToGeneError(str(self) + " != " +
                                     str(self.details + "\t" + str(count)))

  def __str__(self):
    res = (self.refseq + "_exon_" + str(self.exonNumber) + "_0_" +
           str(self.chrom) + "_" + str(self.loc))
    if self.strand == "+":
      res += "_f"
    else:
      res += "_r"
    res += ("\t" + str(self.count))
    return res

  def isIn(self, bedElement):
    if bedElement.chrom != self.chrom:
      return False
    if bedElement.strand != self.strand:
      return False
    if bedElement.name != self.refseq:
      return False
    if self.loc > bedElement.end or self.loc < bedElement.start:
      return False
    return True


###############################################################################
#           USER INTERFACE, COMMAND-LINE PROCESSING AND DISPATCH              #
###############################################################################

def getUI(prog_name, args):
  programName = prog_name
  longDescription = "Given a tab delineated file with two coloumns, one " +\
                    "for exon name and the other for read count, and a BED " +\
                    "file indicating the locations of genes, produce gene " +\
                    "read counts"
  shortDescription = longDescription

  ui = CLI(programName, shortDescription, longDescription)
  ui.minArgs = 0
  ui.maxArgs = 1
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output resultant RPKM values to this file",
                      required=False, type=str))
  ui.addOption(Option(short="v", long="verbose",
                      description="output additional messages to stderr " +
                                  "about run (default: " +
                                  str(DEFAULT_VERBOSITY) + ")",
                      default=DEFAULT_VERBOSITY, required=False))
  ui.addOption(Option(short="r", long="reference", argName="filename",
                      description="reference for transcripts (not exons)",
                      required=True, type=str))
  ui.addOption(Option(short="h", long="help",
                      description="show this help message ", special=True))
  ui.addOption(Option(short="u", long="test",
                      description="run unit tests ", special=True))

  ui.parseCommandLine(args)
  return ui


def _main(args, prog_name):
  # get options and arguments
  ui = getUI(prog_name, args)

  if ui.optionIsSet("test"):
    # just run unit tests
    unittest.main(argv=[sys.argv[0]])
  elif ui.optionIsSet("help"):
    # just show help
    ui.usage()
  else:
    verbose = (ui.optionIsSet("verbose") is True) or DEFAULT_VERBOSITY

    # get the reference for the transcripts/genes
    reffh = open(ui.getValue("reference"))

    # get output filename/handle
    outFh = sys.stdout
    if ui.optionIsSet("output"):
      outFh = open(ui.getValue("output"), "w")

    # figure out where our input is coming from... and read the counts
    if ui.hasArgument(0):
      inFh = open(ui.getArgument(0))
    else:
      inFh = sys.stdin

    process(inFh, reffh, outFh, verbose)


###############################################################################
#                              MAIN PROGRAM LOGIC                             #
###############################################################################

def makeExonDictionary(exons):
  byChrom = {}
  for e in exons:
    if e.chrom not in byChrom:
      byChrom[e.chrom] = {}
    if e.strand not in byChrom[e.chrom]:
      byChrom[e.chrom][e.strand] = {}
    if e.refseq not in byChrom[e.chrom][e.strand]:
      byChrom[e.chrom][e.strand][e.refseq] = []
    byChrom[e.chrom][e.strand][e.refseq].append(e)
  return byChrom


def readExonCounts(fh):
  if type(fh).__name__ == "str":
    fh = open(fh)
  res = []
  for l in fh:
    l = l.strip()
    if l == "":
      continue
    parts = l.split("\t")
    if len(parts) == 2:
      # guess this is a simple two-col file with exon name and count
      details, count = parts
    elif len(parts) == 6:
      # guess this is a BED file
      details, count = parts[3], parts[4]
    else:
      raise readCountExonToGeneError("unknown file format for exon " +
                                     "read counts")
    res.append(Exon(details, count))
  return res


def process(infh, reffh, outfh, verbose=False, debug=False):
  # first read the exons and produce a hash so we can look them up quickly..
  exonHash = makeExonDictionary(readExonCounts(infh))

  # just some debugging messages.. make sure we built the hash properly.
  if debug:
    sys.stderr.write("exonHash is: \n")
    for chrom in exonHash:
      sys.stderr.write(chrom + "\n")
      for strand in exonHash[chrom]:
        sys.stderr.write("\t" + strand + "\n")
        for refseq in exonHash[chrom][strand]:
          sys.stderr.write("\t\t" + refseq + "\n")
          for exon in exonHash[chrom][strand][refseq]:
            sys.stderr.write("\t\t\t" + str(exon) + "\n")

  # now process the genes/transcripts
  for t in BEDIterator(reffh):
    if debug:
      sys.stderr.write("processing " + str(t) + "\n")
    hits = []
    count = 0
    try:
      hits = exonHash[t.chrom][t.strand][t.name]
    except KeyError:
      pass
    if debug:
      sys.stderr.write("hits are: \t" +
                       "\n\t\t".join([str(s) for s in hits]) + "\n")
    for h in hits:
      if h.isIn(t):
        count += h.count
    t.score = count
    outfh.write(str(t) + "\n")


###############################################################################
#                               UNIT TESTS                                    #
###############################################################################

class ExonToGeneReadCountTests(unittest.TestCase):
  """
    Unit tests for converting exon read counts to gene read counts
  """

  def setUp(self):
    self.genes =\
        "\t".join(["chr1", "10", "80", "NM_001", "0", "+"]) + "\n" +\
        "\t".join(["chr1", "20", "120", "NM_002", "0", "-"]) + "\n" +\
        "\t".join(["chr1", "45", "60", "NM_002", "0", "+"]) + "\n" +\
        "\t".join(["chr1", "90", "125", "NM_001", "0", "+"]) + "\n" +\
        "\t".join(["chr2", "20", "75", "NM_002", "0", "+"]) + "\n" +\
        "\t".join(["chr2", "45", "85", "NM_003", "0", "-"]) + "\n"

    self.exonCounts =\
        "NM_001_exon_1_0_chr1_10_f" + "\t" + "10" + "\n" +\
        "NM_001_exon_2_0_chr1_40_f" + "\t" + "10" + "\n" +\
        "NM_001_exon_3_0_chr1_70_f" + "\t" + "10" + "\n" +\
        "NM_002_exon_1_0_chr1_20_r" + "\t" + "20" + "\n" +\
        "NM_002_exon_2_0_chr1_45_r" + "\t" + "20" + "\n" +\
        "NM_002_exon_3_0_chr1_100_r" + "\t" + "20" + "\n" +\
        "NM_002_exon_1_0_chr1_45_f" + "\t" + "30" + "\n" +\
        "NM_001_exon_1_0_chr1_90_f" + "\t" + "40" + "\n" +\
        "NM_001_exon_2_0_chr1_120_f" + "\t" + "40" + "\n" +\
        "NM_002_exon_1_0_chr2_20_f" + "\t" + "10" + "\n" +\
        "NM_002_exon_2_0_chr2_60_f" + "\t" + "10" + "\n" +\
        "NM_003_exon_1_0_chr2_45_r" + "\t" + "20" + "\n" +\
        "NM_003_exon_2_0_chr2_77_r" + "\t" + "20" + "\n"

    self.expectedAns =\
        ["\t".join(["chr1", "10", "80", "NM_001", "30", "+"]) + "\n",
         "\t".join(["chr1", "20", "120", "NM_002", "60", "-"]) + "\n",
         "\t".join(["chr1", "45", "60", "NM_002", "30", "+"]) + "\n",
         "\t".join(["chr1", "90", "125", "NM_001", "80", "+"]) + "\n",
         "\t".join(["chr2", "20", "75", "NM_002", "20", "+"]) + "\n",
         "\t".join(["chr2", "45", "85", "NM_003", "40", "-"]) + "\n"]

  def testSimpleKeys(self):
    debug = False
    infh = DummyInputStream(self.exonCounts)
    inref = DummyInputStream(self.genes)
    outfh = DummyOutputStream()

    process(infh, inref, outfh, verbose=False, debug=debug)
    gotOutput = outfh.itemsWritten()

    if debug:
      print "expected -------"
      for e in self.expectedAns:
        print e
      print "got ------------"
      for e in gotOutput:
        print e

    assert(self.expectedAns == gotOutput)


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
  try:
    _main(sys.argv[1:], sys.argv[0])
  except Exception as e:
    sys.stderr.write("An unexpected exception occurred. Details: " + str(e) +
                     " \n")
