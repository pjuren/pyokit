#!/usr/bin/python

"""
  Date of Creation: 29th May 2010

  XPIPE Component
  Description:   Given a file A, consisting of mapped reads, and a set of
                 files B, to filter against, remove any read from A that
                 also appears in any file in B.

  Copyright (C) 2010
  University of Southern California,
  Philip J. Uren,
  Jin H. Park,
  Andrew D. Smith

  Authors: Emad Bahrami-Samani, Philip J. Uren, Jin H. Park, Andrew D. Smith

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

# python imports
import sys
import unittest
import StringIO

# for testing
import mock

# pyokit imports
from pyokit.interface.cli import CLI, Option
from pyokit.io.bedIterators import BEDUniqueIterator
from pyokit.util.testing import build_mock_open_side_effect


###############################################################################
#                                  CONSTANTS                                  #
###############################################################################

MAX_NEEDED_BED_FIELDS = 6  # need at least this many; we'll drop the rest.
DEFAULT_VERBOSITY = False


###############################################################################
#                               USER INTERFACE                                #
###############################################################################

def getUI(prog_name, args):
  longDescription = "Given two files consisting of mapped reads, remove " +\
                    "the reads that are in both mapped in exonic regions " +\
                    "and junction regions "
  shortDescription = "Filter ambigious reads"

  ui = CLI(prog_name, shortDescription, longDescription)
  ui.minArgs = -1
  ui.maxArgs = -1
  ui.addOption(Option(short="f", long="first", argName="filename",
                      description="BED file with first end reads",
                      required=True, type=str))
  ui.addOption(Option(short="s", long="second", argName="filename",
                      description="BED file with second end reads",
                      required=True, type=str))
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="whitespace separated list of output" +
                                  " filenames, wrap in commas. Expect two" +
                                  " filenames here, one for each input.",
                      required=True, type=str))
  ui.addOption(Option(short="v", long="verbose",
                      description="outputs additional messages to stdout " +
                                  "about run (default: " +
                                  str(DEFAULT_VERBOSITY) + ")",
                      default=DEFAULT_VERBOSITY, required=False))
  ui.addOption(Option(short="b", long="best", description="keep clearly " +
                      "better ambiguously mapped reads", required=False))
  ui.addOption(Option(short="h", long="help",
                      description="show this help message", special=True))
  ui.addOption(Option(short="u", long="test",
                      description="run unit tests", special=True))

  ui.parseCommandLine(args)
  return ui


###############################################################################
#                     COMMAND LINE PROCESSING AND DISPATCH                    #
###############################################################################

def _main(args, prog_name):
  # get options and arguments
  ui = getUI(prog_name, args)
  if ui.optionIsSet("help"):
    # just show help
    ui.usage()
  elif ui.optionIsSet("test"):
    # just run unit tests
    unittest.main(argv=[sys.argv[0]])
  else:
    verbose = (ui.optionIsSet("verbose") is True) or DEFAULT_VERBOSITY

    # get input file handles
    infh1 = open(ui.getValue("first"))
    infh2 = open(ui.getValue("second"))

    # get output file handles
    parts = ui.getValue("output").split()
    if len(parts) != 2:
      sys.stderr.write("expected two output filenames, got" +
                       str(len(parts)))
      sys.exit(1)
    outfh1, outfh2 = open(parts[0], "w"), open(parts[1], "w")

    # only filter worse mappings?
    best = ui.optionIsSet("best")

    ambigFilter(infh1, infh2, outfh1, outfh2, verbose, best)


###############################################################################
#                             MAIN PROGRAM LOGIC                              #
###############################################################################

def ambigFilter(in_fh1, in_fh2, out_fh1, out_fh2, verbose=False, best=False):
  """
    @summary: take reads from in_fh1 and output to out_fh1 if they don't
              appear also in in_fh2 (ditto for in_fh2)

    @param in_fh1: BED formated stream of reads
    @param in_fh2: BED formated stream of reads
    @param out_fh1: Output reads that pass from in_fh1 to this stream
    @param out_fh2: Output reads that pass from in_fh2 to this stream
    @param verbose: output additional messages to sys.stderr if True
    @param best: Given two items that have the same name, try to output the one
                 with the best score

    @return: None (out streams have BED format)
  """
  for r1, r2 in BEDUniqueIterator(in_fh1, in_fh2, verbose, best, dropAfter=6):
    if r1 is not None:
      out_fh1.write(str(r1) + "\n")
    if r2 is not None:
      out_fh2.write(str(r2) + "\n")


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestRemDupsBed(unittest.TestCase):
  """..."""

  def setUp(self):
    self.sfile1 = "chr1 \t 1 \t 2 \t read1\n" +\
                  "chr2 \t 3 \t 4 \t read4\n" +\
                  "chr4 \t 3 \t 3 \t read5\n" +\
                  "chr4 \t 3 \t 3 \t read6"
    self.sfile2 = "chr1 \t 3 \t 4 \t read3\n" +\
                  "chr3 \t 2 \t 4 \t read5"

    self.sfile3 = "chr1 \t 1 \t 2 \t read1 \t 1\n" +\
                  "chr2 \t 3 \t 4 \t read4 \t 0\n" +\
                  "chr4 \t 3 \t 3 \t read5 \t 1\n" +\
                  "chr4 \t 3 \t 3 \t read6 \t 1"
    self.sfile4 = "chr1 \t 3 \t 4 \t read3 \t 2\n" +\
                  "chr3 \t 2 \t 4 \t read5 \t 2"

    self.f_map = {"sfile1": self.sfile1, "sfile2": self.sfile2,
                  "sfile3": self.sfile3, "sfile4": self.sfile4}

  @mock.patch('__builtin__.open')
  def test_basic(self, mock_open):
    outfh1 = StringIO.StringIO()
    outfh2 = StringIO.StringIO()
    streams = {"out1.dat": outfh1, "out2.dat": outfh2}
    mock_open.side_effect = build_mock_open_side_effect(self.f_map, streams)

    _main(["-f", "sfile1", "-s", "sfile2", "-o", "out1.dat out2.dat"],
          sys.argv[0])

    expectedOutput1 = "\n".join(["\t".join(["chr1", "1", "2", "read1"]),
                                 "\t".join(["chr2", "3", "4", "read4"]),
                                 "\t".join(["chr4", "3", "3", "read6"])])
    expectedOutput2 = "\t".join(["chr1", "3", "4", "read3"])

    gotOutput1 = outfh1.getvalue()
    gotOutput2 = outfh2.getvalue()

    self.assertTrue(gotOutput1 == (expectedOutput1 + "\n"))
    self.assertTrue(gotOutput2 == (expectedOutput2 + "\n"))

  @mock.patch('__builtin__.open')
  def test_best_match(self, mock_open):
    outfh1 = StringIO.StringIO()
    outfh2 = StringIO.StringIO()
    streams = {"out1.dat": outfh1, "out2.dat": outfh2}
    mock_open.side_effect = build_mock_open_side_effect(self.f_map, streams)

    _main(["-b", "-f", "sfile3", "-s", "sfile4", "-o", "out1.dat out2.dat"],
          sys.argv[0])

    expectedOutput1 = "\n".join(["\t".join(["chr1", "1", "2", "read1", "1"]),
                                 "\t".join(["chr2", "3", "4", "read4", "0"]),
                                 "\t".join(["chr4", "3", "3", "read5", "1"]),
                                 "\t".join(["chr4", "3", "3", "read6", "1"])])
    expectedOutput2 = "\t".join(["chr1", "3", "4", "read3", "2"])

    gotOutput1 = outfh1.getvalue()
    gotOutput2 = outfh2.getvalue()

    self.assertTrue(gotOutput1 == (expectedOutput1 + "\n"))
    self.assertTrue(gotOutput2 == (expectedOutput2 + "\n"))


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
  try:
    _main(sys.argv[1:], sys.argv[0])
  except Exception as e:
    sys.stderr.write("An unexpected exception occurred. Details: " + str(e) +
                     " \n")
