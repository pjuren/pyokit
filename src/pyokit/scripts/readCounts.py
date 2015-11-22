#! /usr/bin/env python

"""
  Count the number of reads overlapping a set of genomic regions and output the
  counts as BED format -- one line per input region with the count of reads
  as the score field.

  Copyright (C) 2015
  University of Southern California,
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

# stndard python imports
import sys
import unittest

# pyokit imports
from pyokit.datastruct.genomicInterval import bucketIterator
from pyokit.interface.cli import CLI, Option
from pyokit.io.bedIterators import BEDIterator


###############################################################################
#                               USER INTERFACE                                #
###############################################################################

def getUI(prog_name, args):
  longDescription = "..."
  shortDescription = longDescription

  ui = CLI(prog_name, shortDescription, longDescription)
  ui.minArgs = 2
  ui.maxArgs = 2
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output resultant reads to these files",
                      required=False, type=str))
  ui.addOption(Option(short="v", long="verbose",
                      description="output additional messages to stderr " +
                                  "about run.", required=False))
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
    verbose = ui.optionIsSet("verbose")

    # get input handles for reads file and regions file
    read_fh = None
    roi_fh = None
    if ui.hasArgument(0):
      read_fh = open(ui.getArgument(0))
    if ui.hasArgument(1):
      roi_fh = open(ui.getArgument(1))

    # make output handle
    outfh = sys.stdout
    if ui.optionIsSet("output"):
      outfh = open(ui.getValue("output"), "w")

    count(read_fh, roi_fh, outfh, verbose)


###############################################################################
#                             MAIN SCRIPT LOGIC                               #
###############################################################################

def count(reads_fh, rois_fh, outfh, verbose=False):
  for region, reads in bucketIterator(BEDIterator(reads_fh, verbose=verbose),
                                      BEDIterator(rois_fh, verbose=verbose)):
    region.score = len(reads)
    outfh.write(str(region) + "\n")


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class ConvertJunctionsUnitTests(unittest.TestCase):
  def setUp(self):
   pass

###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
  try:
    _main(sys.argv[1:], sys.argv[0])
  except Exception as e:
    msg = str(e).decode('string_escape')
    sys.stderr.write("An unexpected exception occurred. Details: " + str(msg) +
                     " \n")
