"""
Date of Creation: 28st August 2015.

Description:    Join two data tables using a key field in both

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

import sys
import os
import unittest

from interface.cli import CLI, Option
from util.progressIndicator import ProgressIndicator

###############################################################################
#                                 CONSTANTS                                   #
###############################################################################

DEFAULT_VERBOSITY = False


###############################################################################
#                               USER INTERFACE                                #
###############################################################################

def getUI():
  programName = os.path.basename(sys.argv[0])
  longDescription = "join"
  shortDescription = "join"

  ui = CLI(programName, shortDescription, longDescription)
  ui.minArgs = 0
  ui.maxArgs = 1
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output resultant reads to these files",
                      required=False, type=str))
  ui.addOption(Option(short="n", long="names", argName="filename",
                      description="file with reference names and alternative names",
                      required=True, type=str))
  ui.addOption(Option(short="f", long="field", argName="number",
                      description="field in input that contains the reference name",
                      default=1, required=False, type=int))
  ui.addOption(Option(short="r", long="reffield", argName="number",
                      description="field in reference that contains the reference name",
                      default=1, required=False, type=int))
  ui.addOption(Option(short="a", long="add", argName="number",
                      description="field in reference that is to be added",
                      default=1, required=True, type=int))
  ui.addOption(Option(short="e", long="header", argName="fieldHead",
                      description="if set, the first line in the input file is assumed " +\
                                  "to be a header. The new field added will have the " +\
                                  "header value given to this option",
                      required=False, type=str))
  ui.addOption(Option(short="h", long="help",
                      description="show this help message "))
  ui.addOption(Option(short="v", long="verbose",
                      description="outputs additional messages to stdout " +\
                                  "about run (default: " + str(DEFAULT_VERBOSITY) + ")",
                      default=DEFAULT_VERBOSITY, required=False))
  ui.addOption(Option(short="u", long="test",
                      description="run unit tests "))

  ui.parseCommandLine(sys.argv[1:])
  return ui


def _main():
  # get options and arguments
  ui = getUI()

  # just run unit tests
  if ui.optionIsSet("test"):
    unittest.main(argv=[sys.argv[0]])
    sys.exit()

  # just show help
  if ui.optionIsSet("help"):
    ui.usage()
    sys.exit()
  verbose = (ui.optionIsSet("verbose") is True) or DEFAULT_VERBOSITY

  # get input handle
  infh = sys.stdin
  if ui.hasArgument(0):
    infh = open(ui.getArgument(0))

  # make output handle
  outfh = sys.stdout
  if ui.optionIsSet("output"):
    outfh = open(ui.getValue("output"), "w")

  # get reference field
  inputfield = 1
  if ui.optionIsSet("field"):
    inputfield = ui.getValue("field")

  # get reference field
  reffield = 1
  if ui.optionIsSet("reffield"):
    reffield = ui.getValue("reffield")

  # get field number in the reference file to be added
  addField = 1
  if ui.optionIsSet("add"):
    addField = ui.getValue("add")

  # get name lookup table handle
  names = ui.getValue("names")

  # does our input have a header? if so we need to
  # know what to name the new field in the output
  newHeader = None
  if ui.optionIsSet("header") :
    newHeader = ui.getValue("header")

  # do our thing..
  process(infh, outfh, loadNames(names, reffield - 1, addField - 1, verbose), inputfield - 1,
          inputfield - 1, newHeader, verbose)


def fileIterator(filehandle, verbose=False):
  if type(filehandle).__name__ == "str":
    filehandle = open(filehandle)
  if verbose:
    try:
      pind = ProgressIndicator(totalToDo=os.path.getsize(filehandle.name),
                               messagePrefix="completed",
                               messageSuffix="of processing " +
                                             filehandle.name)
    except AttributeError:
      sys.stderr.write("BEDIterator -- warning: " +
                       "unable to show progress for stream")
      verbose = False

  for line in filehandle:
    line = line.strip()
    if verbose:
      pind.done = filehandle.tell()
      pind.showProgress()
    if line == "":
      continue
    yield line


def loadNames(fn, reffield, addField, verbose=False):
  names = {}
  for line in fileIterator(fn, verbose=verbose):
    if len(line.split("\t")) < 2:
      continue
    parts = line.split("\t")
    refname = parts[reffield]
    altname = parts[addField]

    if refname not in names:
      names[refname] = []
    names[refname].append(altname)
  return names


def process(infh, outfh, names, inputfield, outputField=None,
            newHeader=None, verbose=False):
  if outputField is None:
    outputField = inputfield

  first = True
  for line in fileIterator(infh, verbose=verbose):
    parts = line.split("\t")

    # process first line as special case?
    if first:
      first = False
      if newHeader is not None:
        parts.insert(outputField, newHeader)
        outfh.write("\t".join(parts) + "\n")
        continue

    refseq = parts[inputfield]
    namestring = ""
    firsti = True
    for r in refseq.split(";"):
      if r in names:
        if not firsti:
          namestring += ","
        firsti = False
        namestring += ",".join(names[r])
    if namestring == "":
      namestring = "UNKNOWN"

    parts.insert(outputField, namestring)
    outfh.write("\t".join(parts) + "\n")

if __name__ == "__main__":
  _main()
