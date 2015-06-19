"""
Date of Creation: 2010.

Description:    Takes a file with a list of p-values and applies Benjamini
                and Hochberg FDR to convert to q-values.

Copyright (C) 2010
University of Southern California,
Philip J. Uren,
Andrew D. Smith

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
import os
import unittest

# rpy2 imports
from rpy2.robjects import r

# pyokit imports
from pyokit.interface.cli import CLI, Option


###############################################################################
#                                 CONSTANTS                                   #
###############################################################################

DEFAULT_VERBOSITY = False
DEFAULT_DELIM = "\t"


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
  ui.minArgs = 0
  ui.maxArgs = 1
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output to given file, else stdout",
                      required=False, type=str))
  ui.addOption(Option(short="f", long="field", argName="number",
                      description="the field that contains the p-values " +
                                  "(indexed from 1, as with unix 'cut' " +
                                  "command)", default=4,
                      required=False, type=int))
  ui.addOption(Option(short="e", long="header",
                      description="take the first line as a header ",
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


################################################################################
#             PRELIMINARY COMMAND LINE PROCESSING AND DISPATCH                 #
################################################################################

def main(args):
  """
  main entry point for the FDR script.

  :param args: the arguments for this script, as a list of string. Should
               already have had things like the script name stripped. That
               is, if there are no args provided, this should be an empty
               list.
  """
  # get options and arguments
  ui = getUI(args)

  # just run unit tests
  if ui.optionIsSet("test"):
    unittest.main(argv=[sys.argv[0]])
    sys.exit()

  # just show help
  if ui.optionIsSet("help"):
    ui.usage()
    sys.exit()
  verbose = (ui.optionIsSet("verbose") == True) or DEFAULT_VERBOSITY

  # header?
  header = ui.optionIsSet("header")

  # get field value
  field = ui.getValue("field") - 1

  # get output handle
  out_fh = sys.stdout
  if ui.optionIsSet("output"):
    out_fh = open(ui.getValue("output"), "w")

  # get input file-handle
  in_fh = sys.stdin
  if ui.hasArgument(0):
    in_fh = open(ui.getArgument(0))

  delim = DEFAULT_DELIM

  # get input...
  if verbose:
    sys.stderr.write("getting input...\n")
  lines = []
  pvals = []
  headLine = None
  for line in in_fh:
    line = line.strip()
    if line == "":
      continue
    if header and headLine is None:
      headLine = line
      continue
    parts = line.split(delim)
    lines.append(parts)
    pvals.append(parts[field])

  # do conversion....
  if verbose:
    sys.stderr.write("converting...\n")
  strs = []
  for i in pvals:
    strs.append(str(i))
  pvalstr = ",".join(strs)
  r("pvals = c(" + pvalstr + ")")
  qvals = r("p.adjust(pvals, method=\"BH\")")

  # output result...
  if verbose:
    sys.stderr.write("outputing...\n")
  if headLine is not None:
    out_fh.write(headLine + "\n")
  for linenum in range(0, len(lines)):
    line = lines[linenum]
    ln = []
    for i in range(0, len(line)):
      if i == field:
        ln.append(str(qvals[linenum]))
      else:
        ln.append(line[i])
    out_fh.write(delim.join(ln) + "\n")

  in_fh.close()
  out_fh.close()


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == "__main__":
    main()
