"""
Date of Creation: 1st July 2015

Description:    Build an index for a file

Copyright (C) 2010
University of Southern California,
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
import sys
import unittest

# pyokit imports
from pyokit.interface.cli import CLI, Option

# support for enumerations in python 2
from enum import Enum


###############################################################################
#                        CONSTANTS AND ENUMERATIONS                           #
###############################################################################

DEFAULT_VERBOSITY = False


class SupportedFileTypes(Enum):

  """An enumeration of the possible file types supported by the script."""

  maf = 1
  reapeat_masker = 2


###############################################################################
#                              USER INTERFACE                                 #
###############################################################################

def getUI(prog_name, args):
  """
  build and return a UI object for this script.

  :param args: raw arguments to parse
  """
  programName = prog_name
  long_description = "Build an index for one or more files."
  short_description = long_description

  ui = CLI(programName, short_description, long_description)
  ui.minArgs = 0
  ui.maxArgs = 1
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output to given file, else stdout",
                      required=False, type=str))
  ui.addOption(Option(short="t", long="type", argName="filename",
                      description="the type of the file. If missing, " +
                      "the script will try to guess the file type. " +
                      "Supported file types are: " +
                      ", ".join([f.name for f in SupportedFileTypes]),
                      default=4, required=False, type=int))
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


################################################################################
#                     COMMAND LINE PROCESSING AND DISPATCH                     #
################################################################################

def main(args=[], prog_name=sys.argv[0]):
  """
  main entry point for the index script.

  :param args:      the arguments for this script, as a list of string. Should
                    already have had the script name stripped. That is, if
                    there are no args provided, this should be an empty list.
  :param prog_name: the name of the script; taken from command line args by
                    default, but you can change it if you want.
  """
  # get options and arguments
  ui = getUI(prog_name, args)

  # just run unit tests
  if ui.optionIsSet("test"):
    unittest.main(argv=[sys.argv[0]])
    sys.exit()

  # just show help
  if ui.optionIsSet("help"):
    ui.usage()
    sys.exit()
  verbose = (ui.optionIsSet("verbose") == True) or DEFAULT_VERBOSITY

  # get input file-handle(s); fall back to stdin if none found.
  in_fhs = [open(x) for x in ui.getArguments()]
  if in_fhs == []:
    if sys.stdin.isatty():
      sys.stderr.write("[NO INPUT FILE FOUND; WAITING FOR INPUT FROM STDIN]")
    in_fhs = [sys.stdin]

  # Get output handle(s). If only one input file, and no output file, output
  # to standard out.
  out_fhs = [sys.stdout]
  if ui.optionIsSet("output"):
    out_fhs = [open(x, "w") for x in ui.getValue("output").split()]
  if len(in_fhs) != len(out_fhs):
    sys.stderr.write("mismatch between number of input files and output files")
    sys.exit(0)


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == "__main__":
    main(sys.argv[1:], sys.argv[0])
