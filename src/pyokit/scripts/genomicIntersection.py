#!/usr/bin/python

"""
  Date of Creation: Dec 2015
  Description:      Get the intersection of genomic regions in one or more
                    files

  Copyright (C) 2010-2015
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

# standard python imports
import os
import sys
import unittest
import StringIO

# for unit testing
import mock

# pyokit imports
from pyokit.interface.cli import CLI, Option
from pyokit.io.bedIterators import BEDIterator
from pyokit.util.testing import build_mock_open_side_effect
from pyokit.datastruct.genomicInterval import regionsIntersection


################################################################################
#                               USER INTERFACE                                 #
################################################################################

def getUI(args):
  """
  build and return a UI object for this script.

  :param args: raw arguments to parse
  """
  programName = os.path.basename(sys.argv[0])
  longDescription = "takes a file with a list of p-values and applies " +\
                    "Benjamini and Hochberg FDR to convert to q-values "
  shortDescription = "takes a file with a list of p-values and applies " +\
                     "Benjamini and Hochberg FDR to convert to q-values "

  ui = CLI(programName, shortDescription, longDescription)
  ui.minArgs = 2
  ui.maxArgs = 2
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output to given file, else stdout",
                      required=False, type=str))
  ui.addOption(Option(short="s", long="stranded",
                      description="treat regions on separate strands as " +
                                  "disjoint, even if they overlap",
                                  required=False))
  ui.addOption(Option(short="v", long="verbose",
                      description="output additional messages to stderr " +
                                  "about run", required=False))
  ui.addOption(Option(short="h", long="help",
                      description="show this help message ", special=True))
  ui.addOption(Option(short="u", long="test",
                      description="run unit tests ", special=True))

  ui.parseCommandLine(args)
  return ui


################################################################################
#                     COMMAND LINE PROCESSING AND DISPATCH                     #
################################################################################

def main(args):
  """
  main entry point for the GenomicIntIntersection script.

  :param args: the arguments for this script, as a list of string. Should
               already have had things like the script name stripped. That
               is, if there are no args provided, this should be an empty
               list.
  """
  # get options and arguments
  ui = getUI(args)

  if ui.optionIsSet("test"):
    # just run unit tests
    unittest.main(argv=[sys.argv[0]])
  elif ui.optionIsSet("help"):
    # just show help
    ui.usage()
  else:
    verbose = ui.optionIsSet("verbose")

    # stranded?
    stranded = ui.optionIsSet("stranded")
    if stranded:
      sys.stderr.write("Sorry, stranded mode hasn't been implemented yet.")
      sys.exit()

    # get output handle
    out_fh = sys.stdout
    if ui.optionIsSet("output"):
      out_fh = open(ui.getValue("output"), "w")

    # get input file-handles -- we know we'll get exactly two, since we
    # specified it in the UI definition
    regions_1 = [x for x in BEDIterator(ui.getArgument(0), verbose=verbose)]
    regions_2 = [x for x in BEDIterator(ui.getArgument(1), verbose=verbose)]

    for r in regionsIntersection(regions_1, regions_2):
      out_fh.write(str(r) + "\n")


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################

class TestGenomicIntIntersection(unittest.TestCase):

  """Unit tests for the GenomicIntIntersection program/script."""

  def setUp(self):
    f1_conts = "\n".join(["\t".join(["chr1", "10", "20", "F1_R1", "0", "+"]),
                          "\t".join(["chr1", "70", "80", "F2_R2", "0", "+"]),
                          "\t".join(["chr2", "05", "10", "F2_R2", "0", "+"]),
                          "\t".join(["chr2", "70", "75", "F2_R2", "0", "+"]),
                          "\t".join(["chr2", "90", "95", "F2_R2", "0", "+"])])
    f2_conts = "\n".join(["\t".join(["chr1", "07", "12", "F1_R1", "0", "+"]),
                          "\t".join(["chr1", "67", "70", "F2_R2", "0", "+"]),
                          "\t".join(["chr1", "75", "85", "F2_R2", "0", "+"]),
                          "\t".join(["chr2", "20", "30", "F2_R2", "0", "+"]),
                          "\t".join(["chr2", "73", "92", "F2_R2", "0", "+"])])
    self.str_d = {"file1.bed": f1_conts, "file2.bed": f2_conts}

  @mock.patch('__builtin__.open')
  def testTwoFiles(self, mock_open):
    outfh = StringIO.StringIO()
    stream_d = {"out.bed": outfh}
    mock_open.side_effect = build_mock_open_side_effect(self.str_d, stream_d)
    main(["-o", "out.bed", "file1.bed", "file2.bed"])
    expect = "\n".join(["\t".join(["chr1", "10", "12", "X", "0", "+"]),
                        "\t".join(["chr1", "75", "80", "X", "0", "+"]),
                        "\t".join(["chr2", "73", "75", "X", "0", "+"]),
                        "\t".join(["chr2", "90", "92", "X", "0", "+"])]) + "\n"
    self.assertEqual(outfh.getvalue(), expect)


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == "__main__":
    main(sys.argv[1:])
