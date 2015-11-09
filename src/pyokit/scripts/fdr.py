"""
Date of Creation: 2010.

Description:    Takes a file with a list of p-values and applies Benjamini
                and Hochberg FDR to convert to q-values.

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
import os
import sys
import mock
import unittest
import StringIO

# pyokit imports
from pyokit.interface.cli import CLI, Option
from pyokit.statistics.multipleHypothesisTesting import correct_pvals


###############################################################################
#                                 CONSTANTS                                   #
###############################################################################

DEFAULT_VERBOSITY = False
DEFAULT_DELIM = "\t"


################################################################################
#                                DATA TABLE CLASS                              #
################################################################################

class DataTable(object):

  """Class representing a simple column-major data table."""

  def __init__(self):
    """Constructor for DataTable objects. Class docstring has param. details."""
    self.header = None
    self.frame = []

  def clear(self):
    """Clear this data frame of all data."""
    self.header = None
    self.frame = []

  def load(self, in_fh, header=False, delimit=None, verbose=False):
    """
    Load this data_table from a stream or file.

    Blank lines in the file are skipped. Any existing values in this dataTable
    object are cleared before loading the new ones.

    :param in_fh:    load from this stream. Can also be a string, in which case
                     we treat it as a filename and attempt to load from that
                     file.
    :param header:   if True, the first row is considered a header
    :param delimit:  delimiter for splitting columns; set to None (default) to
                     split around any whitespace.
    :param verbose:  if True, output progress messages to stderr.
    """
    self.clear()
    if verbose:
      sys.stderr.write("getting input...\n")

    # figure out whether we need to open a file or not
    in_strm = in_fh
    if isinstance(in_strm, basestring):
      in_strm = open(in_strm)

    for line in in_strm:
      line = line.strip()
      if line == "":
        continue
      if header and self.header is None:
        self.header = line.split(delimit)
        continue
      parts = line.split(delimit)
      if self.frame != [] and len(parts) != len(self.frame):
        raise IOError("Cannot handle ragged data frames")
      while len(self.frame) < len(parts):
        self.frame.append([])
      for i in range(0, len(parts)):
        self.frame[i].append(parts[i])

  def write(self, strm, delim, verbose=False):
    """
    Write this data frame to a stream or file.

    :param strm:     stream to write to; can also be a string, in which case we
                     treat it as a filename and to open that file for writing
                     to.
    :param delim:    delimiter to use between columns.
    :param verbose:  if True, output progress messages to stderr.
    """
    if verbose:
      sys.stderr.write("outputing...\n")

    # figure out whether we need to open a file or not
    out_strm = strm
    if isinstance(out_strm, basestring):
      out_strm = open(out_strm)

    if self.header is not None:
      out_strm.write(delim.join(self.header))
    max_col_len = len(max(self.frame, key=len))
    for i in range(0, max_col_len):
      for j in range(0, len(self.frame)):
        if j != 0:
          out_strm.write(delim)
        out_strm.write(str(self.frame[j][i]))
      out_strm.write("\n")


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
#                     COMMAND LINE PROCESSING AND DISPATCH                     #
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

  if ui.optionIsSet("test"):
    # just run unit tests
    unittest.main(argv=[sys.argv[0]])
  elif ui.optionIsSet("help"):
    # just show help
    ui.usage()
  else:
    verbose = (ui.optionIsSet("verbose") is True) or DEFAULT_VERBOSITY

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

    # load data, do conversion, write out results.
    data_table = DataTable()
    data_table.load(in_fh, header, delim, verbose)
    data_table.frame[field] =\
        correct_pvals(data_table.frame[field], verbose=verbose)
    data_table.write(out_fh, delim, verbose)


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################

class TestFDR(unittest.TestCase):

  """Unit tests for the FDR program/script."""

  def setUp(self):
    """Set up some data for use in the tests."""
    # uniformly distributed numbers between 0 and 1, like p-values...
    self.unif_n_1 = "0.2547 \n 0.6625 \n 0.0331 \n 0.3938 \n 0.0822 \n" +\
                    "0.6728 \n 0.8484 \n 0.9453 \n 0.1612 \n 0.1747 \n" +\
                    "0.8829 \n 0.1905 \n 0.1736 \n 0.6905 \n 0.3057 \n" +\
                    "0.4780 \n 0.1908 \n 0.8138 \n 0.3669 \n 0.9748 \n" +\
                    "0.4855 \n 0.3297 \n 0.9502 \n 0.4268 \n 0.2953 \n" +\
                    "0.1949 \n 0.0690 \n 0.4817 \n 0.4864 \n 0.5723 \n" +\
                    "0.2864 \n 0.8836 \n 0.5056 \n 0.6863 \n 0.4898 \n" +\
                    "0.0140 \n 0.5309 \n 0.0010 \n 0.8609 \n 0.3165 \n" +\
                    "0.7512 \n 0.3589 \n 0.2014 \n 0.7035 \n 0.8815 \n" +\
                    "0.8289 \n 0.2166 \n 0.4445 \n"

    # self.unif_n_1 corrected using BH method.
    self.unif_n_1_corr = [0.8089600, 0.9126486, 0.5296000, 0.8089600,
                          0.7891200, 0.9126486, 0.9425067, 0.9704170,
                          0.7997538, 0.7997538, 0.9425067, 0.7997538,
                          0.7997538, 0.9126486, 0.8089600, 0.8089600,
                          0.7997538, 0.9425067, 0.8089600, 0.9748000,
                          0.8089600, 0.8089600, 0.9704170, 0.8089600,
                          0.8089600, 0.7997538, 0.7891200, 0.8089600,
                          0.8089600, 0.8584500, 0.8089600, 0.9425067,
                          0.8089600, 0.9126486, 0.8089600, 0.3360000,
                          0.8220387, 0.0480000, 0.9425067, 0.8089600,
                          0.9425067, 0.8089600, 0.7997538, 0.9126486,
                          0.9425067, 0.9425067, 0.7997538, 0.8089600]

  @mock.patch('__builtin__.open')
  def test_no_header_simple(self, mock_open):
    """
    A simple test of basic functionality whith a small input "file".

    The file has no header, it has only one column and that column contains
    the p-values to be adjusted.
    """
    dummy_output = StringIO.StringIO()

    # set up our mock for open
    def open_side_effect(*args, **kwargs):
      if args[0] == "input.txt":
        return StringIO.StringIO(self.unif_n_1)
      elif args[0] == "output.txt":
        return dummy_output
      raise IOError("No such file")
    mock_open.side_effect = open_side_effect

    main(["-f", "1", "-o", "output.txt", "input.txt"])
    res = [float(x) for x in dummy_output.getvalue().split("\n")
           if x.strip() != ""]
    self.assertEqual(len(self.unif_n_1_corr), len(res))
    for i in range(0, len(res)):
      self.assertAlmostEqual(self.unif_n_1_corr[i], res[i])


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == "__main__":
    main(sys.argv[1:])
