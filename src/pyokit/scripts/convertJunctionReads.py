#!/usr/bin/python

"""
  Convert a set of junction reads into exon reads.
  There are several choices for how to convert a junction read -
    see command line options below

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
import random
import unittest
import collections
import StringIO

from pyokit.interface.cli import CLI, Option
from pyokit.io.bedIterators import BEDIterator
from pyokit.datastruct.genomicInterval import GenomicInterval, parseBEDString
from pyokit.testing.dummyfiles import DummyInputStream, DummyOutputStream


###############################################################################
#                                  CONSTANTS                                  #
###############################################################################

FIRST_EXON = "FIRST_EXON"
SECOND_EXON = "SECOND_EXON"
BOTH_EXONS = "BOTH_EXONS"
BIGGEST_EXON = "BIGGEST_EXON"
FIVE_PRIME_END = "FIVE_PRIME_END"
DEFAULT_SCHEME = BOTH_EXONS

DEFAULT_VERBOSITY = False


###############################################################################
#                               USER INTERFACE                                #
###############################################################################

def getUI(prog_name, args):
  longDescription = "Convert junction reads into exonic reads by picking " +\
                    "an exon to assign them to"
  shortDescription = longDescription

  ui = CLI(prog_name, shortDescription, longDescription)
  ui.minArgs = 0
  ui.maxArgs = 1
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output resultant reads to these files",
                      required=False, type=str))
  ui.addOption(Option(short="v", long="verbose",
                      description="output additional messages to stderr " +
                                  "about run (default: " +
                                  str(DEFAULT_VERBOSITY) + ")",
                      default=DEFAULT_VERBOSITY, required=False))
  ui.addOption(Option(short="s", long="scheme", argName="TYPE",
                      description="scheme for assigning reads to " +
                                  "junctions, one of FIRST_EXON, " +
                                  "SECOND_EXON, BOTH_EXONS, BIGGEST_EXON",
                      default=DEFAULT_SCHEME, required=False, type=str))
  ui.addOption(Option(short="h", long="help",
                      description="show this help message "))
  ui.addOption(Option(short="u", long="test",
                      description="run unit tests "))

  ui.parseCommandLine(args)
  return ui


###############################################################################
#                     COMMAND LINE PROCESSING AND DISPATCH                    #
###############################################################################

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

    # get input handle
    infh = sys.stdin
    if ui.hasArgument(0):
      infh = open(ui.getArgument(0))

    # make output handle
    outfh = sys.stdout
    if ui.optionIsSet("output"):
      outfh = open(ui.getValue("output"), "w")

    # get scheme
    scheme = DEFAULT_SCHEME
    if ui.optionIsSet("scheme"):
      scheme = ui.getValue("scheme")

    processBED(infh, outfh, scheme, verbose)


###############################################################################
#                             MAIN PROGRAM LOGIC                              #
###############################################################################

def processBED(infh, outhandle, scheme, verbose=False):
  for read in BEDIterator(infh, verbose=verbose):
    # split the chrom field to get the genomic indices..
    y = collections.deque(read.chrom.split("_"))
    while len(y) > 5:
      a = y.popleft()
      a += ("_" + y.popleft())
      y.appendleft(a)
    chrom = y[0]
    chrom1SeqStart = int(y[1])
    chrom1SeqEnd = int(y[2])
    chrom2SeqStart = int(y[3])

    # arbitrarily decide the first exon contains the largest portion of
    # the read if both are the same
    firstExon = None
    secondExon = None
    if scheme != SECOND_EXON:
        firstExon = GenomicInterval(chrom, chrom1SeqStart + read.start - 1,
                                    chrom1SeqEnd, read.name, read.score,
                                    read.strand)
    if scheme != FIRST_EXON:
        end = chrom2SeqStart + (read.end - (chrom1SeqEnd - chrom1SeqStart)) - 1
        secondExon = GenomicInterval(chrom, chrom2SeqStart, end, read.name,
                                     read.score, read.strand)

    # we add %1 or %2 to the end of the read names so they can
    # be distinguished later
    if firstExon is not None:
      firstExon.name = firstExon.name + "%1"
    if secondExon is not None:
      secondExon.name = secondExon.name + "%2"

    if (scheme == FIRST_EXON) or \
       (scheme == BIGGEST_EXON and len(firstExon) >= len(secondExon)) or \
       (scheme == FIVE_PRIME_END and read.strand == "+"):
      out = str(firstExon)
    elif (scheme == SECOND_EXON) or \
         (scheme == BIGGEST_EXON and len(secondExon) > len(firstExon)) or \
         (scheme == FIVE_PRIME_END and read.strand == "-"):
      out = str(secondExon)
    elif scheme == BOTH_EXONS:
      out = str(firstExon) + "\n" + str(secondExon)

    # sanity check -- make sure we create a valid output string
    for l in out.split("\n"):
      e = parseBEDString(l)
      if e.chrom.strip() == "":
        raise ValueError(" got an emtpy chrom -> " + str(read))

    # write output
    outhandle.write(out + "\n")


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

def randomString(length, alpha):
  res = ""
  for i in range(0, length):
    idx = int(random.random() * len(alpha))
    res += alpha[idx]
  return res


def randomName(length):
  return randomString(10, "abcdefghijklmnopqrstuvwxyz123456789")


def randomBEDElement(name=None, chrom=None, start=None, end=None, delim="\t",
                     maxIndex=1000000):
  MAX_SCORE = 30

  if name is None:
    name = randomName(10)
  if chrom is None:
    chrom = randomName(10)

  if start is None:
    start = int(random.random() * (maxIndex - 1))
  if end is None:
    end = int(random.random() * (maxIndex - start) + start)
  score = int(random.random() * MAX_SCORE)
  strand = "-"
  if random.random() <= 0.5:
    strand = "+"

  line = delim.join([chrom, str(start), str(end), name, str(score), strand])
  return parseBEDString(line)


def randomJunctionRead(chroms=None, max=sys.maxint):
  if chroms is None:
    chroms = ["chr1", "chr2", "chr3"]

  # pick 4 random chromosomal co-ordinates, ensuring no duplicates
  coords = None
  while coords is None or len(set(coords)) < 4:
    coords = []
    for i in range(0, 4):
      coords.append(random.randint(0, max))
    coords.sort()
    exon1Start = coords[0]
    exon1End = coords[1]
    exon2Start = coords[2]
    exon2End = coords[3]

  lenExon1 = exon1End - exon1Start
  lenExon2 = exon2End - exon2Start

  # pick a chromosome name
  chromName = (random.choice(chroms) + "_" + str(exon1Start) +
               "_" + str(exon1End) + "_" + str(exon2Start) + "_" +
               str(exon2End))

  # pick start and end indices
  start = random.randint(1, lenExon1)
  end = random.randint(lenExon1 + 1, lenExon1 + lenExon2)

  return randomBEDElement(chrom=chromName, start=start, end=end)


class ConvertJunctionsUnitTests(unittest.TestCase):
  def setUp(self):
    self.schemes = [FIRST_EXON, SECOND_EXON, BIGGEST_EXON, BOTH_EXONS]

    NUM_READS = 1000
    MAX_COORD = 100
    reads = [randomJunctionRead(max=MAX_COORD) for i in range(0, NUM_READS)]

    # throw away reads with duplicate names
    filtered = []
    self.names = {}
    for read in reads:
      if read.name in self.names:
        continue
      self.names[read.name] = 1
      filtered.append(read)
    reads = filtered
    self.readLines = [str(e) for e in reads]

    # figure out the answers
    self.lengths = {}
    self.firstChromEnds = {}
    self.firstChromStarts = {}
    self.secondChromStarts = {}
    self.secondChromEnds = {}
    self.readStarts = {}
    self.readEnds = {}
    for read in reads:
      chrom_parts = read.chrom.split("_")
      n = len(chrom_parts)
      self.lengths[read.name] = len(read)
      self.secondChromStarts[read.name] = int(chrom_parts[n - 1])
      self.secondChromStarts[read.name] = int(chrom_parts[n - 2])
      self.firstChromEnds[read.name] = int(chrom_parts[n - 3])
      self.firstChromStarts[read.name] = int(chrom_parts[n - 4])
      self.readStarts[read.name] = read.start
      self.readEnds[read.name] = read.end

  def testInclusion(self):
    """
    if a read appears in the input, it should appear in the output and
    vice-versa. Number of occurances should be the same too (unless we're
    doing BOTH_EXONS, then it should be twice in the output)
    """

    for scheme in self.schemes:
      infh = StringIO.StringIO("\n".join(self.readLines))
      outfh = StringIO.StringIO()
      processBED(infh, outfh, scheme)

      # see what we get..
      outlines = outfh.getvalue().split("\n")
      outlines = [l for l in outlines if l.strip() != ""]
      outnames = [parseBEDString(line).name[:-2] for line in outlines]
      self.assertTrue(set(outnames) == set(self.names.keys()))

      if scheme == BOTH_EXONS:
        len(outnames) / 2 == len(self.names)

  def testBothExonsScheme(self):
    """
    check that output exons have the right size and relative
    indices from chromsome indices
    """

    # run the code...
    infh = StringIO.StringIO("\n".join(self.readLines))
    outfh = StringIO.StringIO()
    processBED(infh, outfh, BOTH_EXONS)

    # see what we get..
    outlines = outfh.getvalue().split("\n")
    outlines = [l for l in outlines if l.strip() != ""]
    for i in range(0, len(outlines), 2):
      first = outlines[i]
      second = outlines[i + 1]

      e1 = parseBEDString(first)
      e2 = parseBEDString(second)

      self.assertTrue(e1.name[:-2] == e2.name[:-2])
      answer = len(e1) + len(e2)
      self.assertTrue(self.lengths[e1.name[:-2]] == answer)
      self.assertTrue(e1.end == self.firstChromEnds[e1.name[:-2]])
      self.assertTrue(e2.start == self.secondChromStarts[e2.name[:-2]])

  def testFirstExonScheme(self):
    """
    check that output exons have the right size and relative indices
    from chromosome indices
    """

    # run the code...
    infh = DummyInputStream(self.readLines)
    outfh = DummyOutputStream()
    processBED(infh, outfh, FIRST_EXON)

    # see what we get..
    outlines = [l.strip() for l in outfh.itemsWritten() if l.strip() != ""]

    for i in range(0, len(outlines)):
      out = outlines[i]
      e1 = parseBEDString(out)

      gotAnswer = len(e1)

      read_start_global = (self.firstChromStarts[e1.name[:-2]] +
                           self.readStarts[e1.name[:-2]])
      expectedAns = self.firstChromEnds[e1.name[:-2]] - read_start_global + 1
      self.assertTrue(gotAnswer == expectedAns)
      self.assertTrue(e1.end == self.firstChromEnds[e1.name[:-2]])

  def testSecondExonScheme(self):
    """
    check that output exons have the right size and relative indices
    from chromosome indices.
    """

    # run the code...
    infh = DummyInputStream(self.readLines)
    outfh = DummyOutputStream()
    processBED(infh, outfh, SECOND_EXON)

    # see what we get..
    outlines = [l.strip() for l in outfh.itemsWritten() if l.strip() != ""]

    for i in range(0, len(outlines)):
      out = outlines[i]
      e2 = parseBEDString(out)

      gotAnswer = len(e2)
      r_len = self.readEnds[e2.name[:-2]] - self.readStarts[e2.name[:-2]]
      glob_s = (self.firstChromStarts[e2.name[:-2]] +
                self.readStarts[e2.name[:-2]])
      expectedAns = r_len - (self.firstChromEnds[e2.name[:-2]] - (glob_s)) - 1
      self.assertTrue(gotAnswer == expectedAns)
      self.assertTrue(e2.start == self.secondChromStarts[e2.name[:-2]])


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
  try:
    _main(sys.argv[1:], sys.argv[0])
  except Exception as e:
    sys.stderr.write("An unexpected exception occurred. Details: " + str(e) +
                     " \n")
