"""
Date of Creation: 2015.

Description:    Collapse the genomic regions in an input stream

Copyright (C) 2015
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
import mock
import unittest
import StringIO

# pyokit imports
from pyokit.datastruct.genomicInterval import collapseRegions
from pyokit.io.bedIterators import BEDIterator
from pyokit.interface.cli import CLI, Option


###############################################################################
#                                 CONSTANTS                                   #
###############################################################################

DEFAULT_VERBOSITY = False
NAME_DELIM = ";"
KEY_DELIM = "\t"


###############################################################################
#                               USER INTERFACE                                #
###############################################################################

def getUI(args, prog_name):
  """
  build and return a UI object for this script.

  :param args: raw arguments to parse
  """
  longDescription = "Collapse the regions in a file of genomic intervals "
  shortDescription = longDescription

  ui = CLI(prog_name, shortDescription, longDescription)
  ui.minArgs = 1
  ui.maxArgs = 1
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output to given file, else stdout",
                      required=False, type=str))
  ui.addOption(Option(short="s", long="stranded",
                      description="collapse positive and negative " +
                                  "strands independently (otherwise strand " +
                                  "is ignored and all output regions are on " +
                                  "positive strand", default=False,
                      required=False))
  ui.addOption(Option(short="a", long="accumulate_names",
                      description="accumulate names of the regions, " +
                                  "separated by semi-colon (otherwise " +
                                  "all output region names will be 'X').",
                      default=False, required=False))
  ui.addOption(Option(short="e", long="exact",
                      description="collapse only regions with exactly the " +
                                  " same genomic location (otherwise any " +
                                  "overlapping regions are collapsed)",
                      default=False, required=False))
  ui.addOption(Option(short="v", long="verbose",
                      description="output additional messages to stderr " +
                                  "about run (default: " +
                                  str(DEFAULT_VERBOSITY) + ")",
                      default=DEFAULT_VERBOSITY, required=False))
  ui.addOption(Option(short="h", long="help",
                      description="show this help message ", special=True))
  ui.addOption(Option(short="u", long="test",
                      description="run unit tests ", special=True))

  ui.parseCommandLine(args)
  return ui


###############################################################################
#                     COMMAND LINE PROCESSING AND DISPATCH                    #
###############################################################################

def main(args, prog_name):
  """
  main entry point for the script.

  :param args: the arguments for this script, as a list of string. Should
               already have had things like the script name stripped. That
               is, if there are no args provided, this should be an empty
               list.
  """
  # get options and arguments
  ui = getUI(args, prog_name)

  # just run unit tests
  if ui.optionIsSet("test"):
    unittest.main(argv=[sys.argv[0]])
    sys.exit()

  # just show help
  if ui.optionIsSet("help"):
    ui.usage()
    sys.exit()
  verbose = (ui.optionIsSet("verbose") is True) or DEFAULT_VERBOSITY

  # how to handle strand, names, and whether to collapse only regions with
  # exactly matching genomic loci?
  stranded = ui.optionIsSet("stranded")
  names = ui.optionIsSet("accumulate_names")
  exact = ui.optionIsSet("exact")

  # get output handle
  out_fh = sys.stdout
  if ui.optionIsSet("output"):
    out_fh = open(ui.getValue("output"), "w")

  # get input file-handle
  in_fh = sys.stdin
  if ui.hasArgument(0):
    in_fh = open(ui.getArgument(0))

  # load data -- TODO at the moment we load everying; need to think about
  # whether it is possible to do this using a single pass of the data, but not
  # loading it all first.
  regions = [x for x in BEDIterator(in_fh, verbose)]

  if exact:
    collapse_exact(regions, stranded, names, out_fh)
  else:
    for x in collapseRegions(regions, stranded, names, verbose):
      out_fh.write(str(x) + "\n")


###############################################################################
#                              COLLAPSING EXACT                               #
###############################################################################

def collapse_exact(regions, stranded, combine_names, out_strm):
  regions_d = {}

  for e in regions:
    key = (KEY_DELIM.join([e.chrom, str(e.start), str(e.end), e.strand]) if
           stranded else KEY_DELIM.join([e.chrom, str(e.start), str(e.end)]))
    if key not in regions_d:
      regions_d[key] = e
      regions_d[key].name = set([regions_d[key].name])
    else:
      if combine_names:
        regions_d[key].name.add(e.name)
  for key in regions_d:
    regions_d[key].name = ";".join(regions_d[key].name)
    out_strm.write(str(regions_d[key]) + "\n")


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################

class TestCollapseRegions(unittest.TestCase):

  """Unit tests for the region collapse script."""

  def setUp(self):
    self.test_in_1 = ["\t".join(["chr1", "10", "20", "R01", "0", "+"]),
                      "\t".join(["chr1", "30", "40", "R02", "1", "+"]),
                      "\t".join(["chr1", "35", "50", "R03", "0", "+"]),
                      "\t".join(["chr1", "45", "65", "R04", "0", "+"]),
                      "\t".join(["chr1", "55", "60", "R05", "3", "-"]),
                      "\t".join(["chr1", "70", "80", "R06", "0", "+"]),
                      "\t".join(["chr1", "75", "95", "R07", "0", "+"]),
                      "\t".join(["chr1", "85", "90", "R08", "1", "-"]),
                      "\t".join(["chr2", "40", "60", "R10", "0", "+"]),
                      "\t".join(["chr3", "10", "20", "R11", "4", "+"]),
                      "\t".join(["chr3", "20", "30", "R12", "0", "-"])]

    self.test_in_2 = ["\t".join(["chr1", "10", "20", "R01", "0", "+"]),
                      "\t".join(["chr1", "10", "20", "R02", "1", "+"]),
                      "\t".join(["chr1", "10", "20", "R01", "0", "+"]),
                      "\t".join(["chr1", "45", "65", "R04", "0", "+"]),
                      "\t".join(["chr1", "55", "60", "R05", "3", "-"]),
                      "\t".join(["chr1", "55", "60", "R06", "0", "+"]),
                      "\t".join(["chr1", "75", "95", "R07", "0", "+"]),
                      "\t".join(["chr1", "85", "90", "R08", "1", "-"]),
                      "\t".join(["chr2", "40", "60", "R10", "0", "+"]),
                      "\t".join(["chr3", "10", "20", "R11", "4", "+"]),
                      "\t".join(["chr3", "10", "20", "R12", "0", "-"])]

  @mock.patch('__builtin__.open')
  def test_maintain_names(self, mock_open):
    """
    Test that names are properly maintinaed when that option is used.
    """

    dummy_output = StringIO.StringIO()

    def open_side_effect(*args, **kwargs):
      if (args[0] == "input.bed"):
        return StringIO.StringIO("\n".join(self.test_in_2))
      elif (args[0] == "output.bed"):
        return dummy_output
      else:
        raise ValueError("Testing failed, unexpected filename: " + args[0])

    mock_open.side_effect = open_side_effect

    # with exact matches, ignore strand
    main(["-a", "-e", "-o", "output.bed", "input.bed"], "regionCollapse")
    r1 = dummy_output.getvalue().strip().split("\n")
    r1.sort()
    expct = ["\t".join(["chr1", "10", "20", "R01;R02", "0", "+"]),
             "\t".join(["chr1", "45", "65", "R04", "0", "+"]),
             "\t".join(["chr1", "55", "60", "R05;R06", "3", "-"]),
             "\t".join(["chr1", "75", "95", "R07", "0", "+"]),
             "\t".join(["chr1", "85", "90", "R08", "1", "-"]),
             "\t".join(["chr2", "40", "60", "R10", "0", "+"]),
             "\t".join(["chr3", "10", "20", "R11;R12", "4", "+"])]
    self.assertEqual(len(r1), len(expct))
    for i in range(0, len(r1)):
      e_lst = expct[i].split("\t")
      r_lst = r1[i].split("\t")
      self.assertEqual(len(e_lst), 6)
      self.assertEqual(len(r_lst), 6)
      for j in [0, 1, 2, 4, 5]:
        self.assertEqual(e_lst[j], r_lst[j])
      self.assertEqual(set(e_lst[3].split(";")), set(r_lst[3].split(";")))

    # with exact matches, separate strands
    # dummy_output = StringIO.StringIO()
    # main(["-a", "-s", "-e", "-o", "output.bed", "input.bed"], "regionCollapse")

    # with overlapping matches, ignore strands
    # dummy_output = StringIO.StringIO()
    # main(["-a", "-o", "output.bed", "input.bed"], "regionCollapse")

    # with overlapping matches, separate strands
    # dummy_output = StringIO.StringIO()
    # main(["-a", "-s", "-o", "output.bed", "input.bed"], "regionCollapse")


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    main(sys.argv[1:], sys.argv[0])
