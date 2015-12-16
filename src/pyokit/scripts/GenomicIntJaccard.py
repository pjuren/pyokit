#!/usr/bin/python

"""
  Date of Creation: 3rd Apr 2010
  Description:      Iterators for processing BED format streams

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

# pyokit imports
from pyokit.interface.cli import CLI, Option
from pyokit.io.bedIterators import BEDIterator
from pyokit.datastruct.genomicInterval import regionsIntersection
from pyokit.datascruct.genomicInterval import collapseRegions


###############################################################################
#                             HELPER FUNCTIONS                                #
###############################################################################
def count_toto_region_size(s):
  """
  sum the size of regions in s; no check for whether they overlap is made.
  """
  tot = 0
  for r in s:
    tot += len(r)


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

    intersection = regionsIntersection(regions_1, regions_2)
    union = collapseRegions(regions_1 + regions_1)

    size_intersection = count_toto_region_size(intersection)
    size_union = count_toto_region_size(union)

    jac = size_intersection / float(size_union)
    print str(jac)


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestGenomicIntJaccard(unittest.TestCase):
  pass


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
