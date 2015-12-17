#!/usr/bin/python

"""
  Date of Creation: Dec 2015
  Description:      Compute the Jaccard index for regions in two files

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

# python imports
import os
import sys
import unittest
import StringIO

# for unit testing
import mock

# pyokit imports
from pyokit.interface.cli import CLI, Option
from pyokit.io.bedIterators import BEDIterator
from pyokit.datastruct.genomicInterval import jaccardIndex
from pyokit.util.testing import build_mock_open_side_effect


###############################################################################
#                              USER INTERFACE                                 #
###############################################################################

def getUI(args):
  """
  build and return a UI object for this script.

  :param args: raw arguments to parse
  """
  programName = os.path.basename(sys.argv[0])
  longDescription = "Compute the Jaccard index for two sets of genomic " +\
                    "intervals."
  shortDescription = longDescription

  ui = CLI(programName, shortDescription, longDescription)
  ui.minArgs = 2
  ui.maxArgs = 2
  ui.addOption(Option(short="v", long="verbose",
                      description="output additional messages to stderr " +
                                  "about run", required=False))
  ui.addOption(Option(short="s", long="stranded",
                      description="treat regions on separate strands as " +
                                  "disjoint, even if they overlap",
                                  required=False))
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
  main entry point for the GenomicIntJaccard script.

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
    stranded = ui.optionIsSet("stranded")

    if stranded:
      sys.stderr.write("Sorry, stranded mode hasn't been implemented yet.")
      sys.exit()

    # we required two input files, so we know these will be present...
    regions_1 = [e for e in BEDIterator(ui.getArgument(0), verbose=verbose)]
    regions_2 = [e for e in BEDIterator(ui.getArgument(1), verbose=verbose)]

    print jaccardIndex(regions_1, regions_2)


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestGenomicIntJaccard(unittest.TestCase):

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
    string_d = {"file1.bed": f1_conts,
                "file2.bed": f2_conts}
    self.mock_open_side_effect = build_mock_open_side_effect(string_d, {})

  @mock.patch('__builtin__.open')
  @mock.patch('sys.stdout', new_callable=StringIO.StringIO)
  def test_simple(self, mock_stdout, mock_open):
    mock_open.side_effect = self.mock_open_side_effect
    main(["file1.bed", "file2.bed"])
    self.assertAlmostEqual(float(mock_stdout.getvalue()), 0.1549296)


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
