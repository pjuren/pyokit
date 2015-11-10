#!/usr/bin/python

"""
  Date of Creation: 3rd Apr 2010
  Description:      Script for generating a profile by piling up a set of
                    regions and counting the number of times another set of
                    intervals overlap these.

  Copyright (C) 2010-2015
  Philip J. Uren,

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
import StringIO
import sys
import unittest

# for testing
import mock

# pyokit imports
from pyokit.interface.cli import CLI
from pyokit.interface.cli import Option
from pyokit.io.bedIterators import intervalTrees
from pyokit.io.bedIterators import BEDIterator
from pyokit.util.testing import build_mock_open_side_effect


###############################################################################
#                                CONSTANTS                                    #
###############################################################################

START = 0
MID = 1
END = 2
WHOLE = 3


###############################################################################
#                             MAIN SCRIPT LOGIC                               #
###############################################################################

def process_anchor_start(regions_fn, to_count_fn, verbose=False):
  """
  :return: list where each element is the number of hits at that location
           relative to the start of the regions.
  """
  res = []
  possible = []
  trees = intervalTrees(to_count_fn, verbose=verbose)
  for region in BEDIterator(regions_fn, verbose=verbose):
    for i in range(0, len(region)):
      while len(possible) <= i:
        possible.append(0)
      possible[i] += 1

    if region.chrom not in trees:
      continue
    hits = trees[region.chrom].intersectingInterval(region.start, region.end)
    for h in hits:
      abs_pos = h.start
      if abs_pos < region.start or abs_pos >= region.end:
        continue
      rel_pos = abs_pos - region.start
      while(len(res) <= rel_pos):
        res.append(0)
      res[rel_pos] += 1

  # normalize
  assert(len(res) <= len(possible))
  for i in range(0, len(res)):
    res[i] = res[i] / float(possible[i])

  return res


###############################################################################
#                                   OUTPUT                                    #
###############################################################################
def write_output(res, out_fh, verbose):
  for i in range(0, len(res)):
    out_fh.write(str(i) + "\t" + str(res[i]) + "\n")


###############################################################################
#                               USER INTERFACE                                #
###############################################################################

def getUI(prog_name, args):
  """Build and return user interface object for this script."""
  longDescription = "Given a set of BED intervals, pile them up so that " +\
                    "their start positions are aligned and then count the " +\
                    "number of other regions that overlap them."
  shortDescription = longDescription

  ui = CLI(prog_name, shortDescription, longDescription)
  # gotta have two args -- MAF dir/file and BED regions.
  # Input by stdin not allowed
  ui.minArgs = 2
  ui.maxArgs = 2
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output to given file, else stdout",
                      required=False, type=str))
  ui.addOption(Option(short="v", long="verbose",
                      description="output additional messages to stderr " +
                                  "about run", required=False))
  ui.addOption(Option(short="h", long="help",
                      description="show this help message ", special=True))
  ui.addOption(Option(short="u", long="test",
                      description="run unit tests ", special=True))

  ui.parseCommandLine(args)
  return ui


###############################################################################
#                     COMMAND LINE PROCESSING AND DISPATCH                    #
###############################################################################

def _main(args, prog_name):
  """Process the command line arguments of this script and dispatch."""
  # get options and arguments
  ui = getUI(prog_name, args)

  if ui.optionIsSet("test"):
    # just run unit tests
    unittest.main(argv=[sys.argv[0]])
  elif ui.optionIsSet("help"):
    # just show help
    ui.usage()
  else:
    verbose = ui.optionIsSet("verbose")

    # know there are exactly two arguments, because we specified this in the
    # UI description.
    regions_fn = ui.getArgument(0)
    to_count_fn = ui.getArgument(1)

    # get output handle
    out_fh = sys.stdout
    if ui.optionIsSet("output"):
      out_fh = open(ui.getValue("output"), "w")

    res = process_anchor_start(regions_fn, to_count_fn, verbose)
    write_output(res, out_fh, verbose)


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestOverlapProfile(unittest.TestCase):

  """Unit tests for this script."""

  def setUp(self):
    self.regions = "chr1 \t 10 \t 20 \t reg1\n" +\
                   "chr1 \t 50 \t 70 \t reg2\n" +\
                   "chr2 \t 15 \t 30 \t reg3\n"
    self.h1 = "chr1 \t 11 \t 14 \t r1\n" +\
              "chr1 \t 15 \t 16 \t r2\n" +\
              "chr1 \t 40 \t 55 \t r3\n" +\
              "chr1 \t 51 \t 53 \t r4\n" +\
              "chr1 \t 55 \t 60 \t r5\n" +\
              "chr1 \t 65 \t 66 \t r6\n" +\
              "chr2 \t 10 \t 11 \t r7\n" +\
              "chr2 \t 16 \t 20 \t r8\n" +\
              "chr2 \t 33 \t 37 \t r9\n"
    self.f_map = {"regions.bed": self.regions, "hits1.bed": self.h1}

  @mock.patch('__builtin__.open')
  def test_anchor_start_norm(self, mock_open):
    outfh = StringIO.StringIO()
    streams = {"out.dat": outfh}
    mock_open.side_effect = build_mock_open_side_effect(self.f_map, streams)

    _main(["-o", "out.dat", "regions.bed", "hits1.bed"], sys.argv[0])
    expect = [[0, 0.0], [1, 1.0], [2, 0.0], [3, 0.0], [4, 0.0], [5, 0.6667],
              [6, 0.0], [7, 0.0], [8, 0.0], [9, 0.0], [10, 0.0], [11, 0.0],
              [12, 0.0], [13, 0.0], [14, 0.0], [15, 1.0]]
    got = outfh.getvalue()
    got = [[float(a) for a in x.split()] for x in got.split("\n")
           if x.strip() != ""]
    self.assertEqual(len(got), len(expect))
    for i in range(0, len(got)):
      self.assertEqual(len(got[i]), 2)
      self.assertAlmostEqual(got[i][0], expect[i][0])
      self.assertAlmostEqual(got[i][1], expect[i][1], places=4)


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    _main(sys.argv[1:], sys.argv[0])
