"""
Date of Creation: 28th August 2015.

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
import mock
import StringIO

from pyokit.interface.cli import CLI, Option
from pyokit.util.progressIndicator import ProgressIndicator

###############################################################################
#                                 CONSTANTS                                   #
###############################################################################

DEFAULT_VERBOSITY = False


###############################################################################
#                            EXCEPTION CLASSES                                #
###############################################################################

class MissingKeyError(Exception):

  """..."""

  def __init__(self, msg):
    """..."""
    self.value = msg

  def __str__(self):
    """Get string representation of this exception."""
    return repr(self.value)


###############################################################################
#                              HELPER FUNCTIONS                               #
###############################################################################

def file_iterator(filehandle, verbose=False):
  """Iterate over a file and yield stripped lines. Optionally show progress."""
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


def load_file(fn_or_strm, key, key_is_field_number, require_unique_key=True,
              delim="\t", verbose=False):
  res = {}
  header = None
  key_field_num = key if key_is_field_number else None
  for line in file_iterator(fn_or_strm, verbose):
    parts = line.split(delim)
    if header is None and not key_is_field_number:
      header = parts
      key_field_num = header.index(key)
      # TODO deal with case where key is not in the header
      # TODO deal with case where key occurs more than once in the header
    else:
      key_val = parts[key_field_num]
      if key_val in res:
        raise ValueError("Oops")
      if key_is_field_number:
        # the entry in the dictionary is a list, minus the key field, in the
        # order they occur.
        res[key_val] = [parts[i] for i in range(0, len(parts))
                        if i != key_field_num]
      else:
        # the entry in the dictionary is another dictionary indexed by
        # the header value
        # TODO proper exception here
        assert(len(parts) == len(header))
        res[key_val] = {}
        for i in range(0, len(parts)):
          if i == key_field_num:
            continue
          else:
            res[key_val][header[i]] = parts[i]
  return res, [x for x in header if x != key] if header is not None else None


###############################################################################
#                               USER INTERFACE                                #
###############################################################################

def getUI(prog_name, args):
  longDescription = "join two data tables using a common field"
  shortDescription = "join two data tables using a common field"

  ui = CLI(prog_name, shortDescription, longDescription)
  ui.minArgs = 2
  ui.maxArgs = 2
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output resultant table to this file. " +
                                  "If omitted, output is to stdout",
                      required=False, type=str))
  ui.addOption(Option(short="a", long="field-one", argName="number/name",
                      description="field in first file to join on. If a " +
                                  "number is given, the file is assumed to " +
                                  "have no header, if a name is given the " +
                                  "file is assumed to have a header. " +
                                  "Columns are indexed from 1. Default is 1",
                      default=1, required=False, type=str))
  ui.addOption(Option(short="b", long="field-two", argName="number/name",
                      description="field in second file to join on. If a " +
                                  "number is given, the file is assumed to " +
                                  "have no header, if a name is given the " +
                                  "file is assumed to have a header. " +
                                  "Columns are indexed from 1. Default is 1",
                      default=1, required=False, type=str))
  ui.addOption(Option(short="m", long="missing", argName="value",
                      description="populate missing fields with this value. " +
                                  "If not provided, missing fields cause " +
                                  "the program to exit with an error.",
                      required=False, type=str))
  ui.addOption(Option(short="h", long="help",
                      description="show this help message ", special=True))
  ui.addOption(Option(short="v", long="verbose",
                      description="outputs additional messages to stdout " +
                                  "about run (default: " +
                                  str(DEFAULT_VERBOSITY) + ")",
                      default=DEFAULT_VERBOSITY, required=False))
  ui.addOption(Option(short="u", long="test",
                      description="run unit tests ", special=True))

  ui.parseCommandLine(args)
  return ui


###############################################################################
#                     COMMAND LINE PROCESSING AND DISPATCH                    #
###############################################################################

def _main(args, prog_name):
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
  verbose = (ui.optionIsSet("verbose") is True) or DEFAULT_VERBOSITY

  # get input handles -- we required two args, so we know these will be here.
  # TODO add unit test for this
  infh1 = open(ui.getArgument(0))
  infh2 = open(ui.getArgument(1))

  # make output handle
  outfh = sys.stdout
  if ui.optionIsSet("output"):
    outfh = open(ui.getValue("output"), "w")

  # get key field in file one
  key_one = 1
  key_one_is_field_number = False
  if ui.optionIsSet("field-one"):
    key_one = ui.getValue("field-one")
    try:
      key_one = int(key_one) - 1
      key_one_is_field_number = True
    except ValueError:
      pass

  # get key field in file two
  key_two = 1
  key_two_is_field_number = False
  if ui.optionIsSet("field-two"):
    key_two = ui.getValue("field-two")
    try:
      key_two = int(key_two) - 1
      key_two_is_field_number = True
    except ValueError:
      pass

  # get missing value
  missing_val = None
  if ui.optionIsSet("missing"):
    missing_val = ui.getValue("missing")

  # do our thing..
  process(infh1, infh2, outfh, key_one, key_one_is_field_number, key_two,
          key_two_is_field_number, missing_val, verbose)


###############################################################################
#                             MAIN PROGRAM LOGIC                              #
###############################################################################

def process(infh1, infh2, outfh, key_one, key_one_is_field_number, key_two,
            key_two_is_field_number, missing_val, verbose=False):
  delim = "\t"

  f2_dictionary, f2_header = load_file(infh2, key_two, key_two_is_field_number,
                                       verbose=verbose)
  f1_header = None
  key_field_num = key_one if key_one_is_field_number else None
  for line in file_iterator(infh1, verbose):
    parts = line.split(delim)
    if f1_header is None and not key_one_is_field_number:
      f1_header = parts
      key_field_num = f1_header.index(key_one)
      # TODO deal with case where key is not in the header
      # TODO deal with case where key occurs more than once in the header

      # TODO exception if header is missing and we have no missing val
      # TODO join missing vals otherwise..
      outfh.write(delim.join(f1_header) + delim + delim.join(f2_header) + "\n")
    else:
      key_val = parts[key_field_num]
      outfh.write(delim.join(parts) + delim)

      if f2_header is not None:
        # f1_dictionary entry is another dictionary indexed by header value
        first = True
        for h_val in f2_header:
          if h_val == key_two:
            continue

          if first:
            first = False
          else:
            outfh.write(delim)
          outfh.write(f2_dictionary[key_val][h_val])
        outfh.write("\n")
      else:
        # f1_dictionary entry is a list
        outfh.write(delim.join(f2_dictionary[key_val]) + "\n")


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestConservationProfileIndvFiles(unittest.TestCase):

  """Unit tests for this script."""

  def setUp(self):
    # matches cleanly with test_file_two
    self.test_headr_one = "\t".join(["AA", "BB", "CC", "DD"])
    self.test_file_one = ["\t".join(["A1", "B1", "C1", "D1"]),
                          "\t".join(["A2", "B2", "C2", "D2"]),
                          "\t".join(["A3", "B3", "C3", "D3"]),
                          "\t".join(["A4", "B4", "C4", "D4"]),
                          "\t".join(["A5", "B5", "C5", "D5"])]

    # matches cleanly with test_file_one
    self.test_headr_two = "\t".join(["XX", "BX", "YY", "ZZ"])
    self.test_file_two = ["\t".join(["X1", "B1", "Y1", "Z1"]),
                          "\t".join(["X2", "B2", "Y2", "Z2"]),
                          "\t".join(["X3", "B3", "Y3", "Z3"]),
                          "\t".join(["X4", "B4", "Y4", "Z4"]),
                          "\t".join(["X5", "B5", "Y5", "Z5"])]

    # macthes with test_file_one, but the keys are not unique in col BB
    self.test_headr_three = "\t".join(["XX", "BX", "YY", "ZZ"])
    self.test_file_three = ["\t".join(["X1", "B1", "Y1", "Z1"]),
                            "\t".join(["X2", "B1", "Y2", "Z2"]),
                            "\t".join(["X3", "B3", "Y3", "Z3"]),
                            "\t".join(["X4", "B4", "Y4", "Z4"]),
                            "\t".join(["X5", "B5", "Y5", "Z5"])]

    # macthes with test_file_two, but the keys are not unique in col BB
    self.test_headr_four = "\t".join(["AA", "BB", "CC", "DD"])
    self.test_file_four = ["\t".join(["A1", "B1", "C1", "D1"]),
                           "\t".join(["A2", "B1", "C2", "D2"]),
                           "\t".join(["A3", "B3", "C3", "D3"]),
                           "\t".join(["A4", "B3", "C4", "D4"]),
                           "\t".join(["A5", "B5", "C5", "D5"])]

    # matches with test_file_two, but some keys in BB are missing
    self.test_headr_five = "\t".join(["AA", "BB", "CC", "DD"])
    self.test_file_five = ["\t".join(["A1", "B1", "C1", "D1"]),
                           "\t".join(["A2", "B2", "C2", "D2"]),
                           "\t".join(["A3", "B3", "C3", "D3"]),
                           "\t".join(["A4", "B4", "C4", "D4"]),
                           "\t".join(["A5", "B5", "C5", "D5"])]

  @mock.patch('__builtin__.open')
  def test_simple_headerless_join(self, mock_open):
    """Test joining two files where all keys are present in both, no header."""
    out_strm = StringIO.StringIO()

    def open_side_effect(*args, **kwargs):
      if args[0] == "one.dat":
        return StringIO.StringIO("\n".join(self.test_file_one))
      if args[0] == "two.dat":
        return StringIO.StringIO("\n".join(self.test_file_two))
      if args[0] == "out.dat":
        return out_strm
      raise IOError("No such file")

    mock_open.side_effect = open_side_effect

    _main(["-a", "2", "-b", "2", "-o", "out.dat", "one.dat", "two.dat"],
          "join")
    expect = ["\t".join(["A1", "B1", "C1", "D1", "X1", "Y1", "Z1"]),
              "\t".join(["A2", "B2", "C2", "D2", "X2", "Y2", "Z2"]),
              "\t".join(["A3", "B3", "C3", "D3", "X3", "Y3", "Z3"]),
              "\t".join(["A4", "B4", "C4", "D4", "X4", "Y4", "Z4"]),
              "\t".join(["A5", "B5", "C5", "D5", "X5", "Y5", "Z5"])]
    self.assertEqual("\n".join(expect) + "\n", out_strm.getvalue())

  @mock.patch('__builtin__.open')
  def test_simple_header_join(self, mock_open):
    """Test joining two files where all keys present in both, with header."""
    out_strm = StringIO.StringIO()

    def open_side_effect(*args, **kwargs):
      if args[0] == "one.dat":
        return StringIO.StringIO(self.test_headr_one + "\n" +
                                 "\n".join(self.test_file_one))
      if args[0] == "two.dat":
        return StringIO.StringIO(self.test_headr_two + "\n" +
                                 "\n".join(self.test_file_two))
      if args[0] == "out.dat":
        return out_strm
      raise IOError("No such file")

    mock_open.side_effect = open_side_effect

    _main(["-a", "BB", "-b", "BX", "-o", "out.dat", "one.dat", "two.dat"],
          "join")
    expect = "\t".join(["AA", "BB", "CC", "DD", "XX", "YY", "ZZ"]) + "\n" +\
             "\n".join(["\t".join(["A1", "B1", "C1", "D1", "X1", "Y1", "Z1"]),
                        "\t".join(["A2", "B2", "C2", "D2", "X2", "Y2", "Z2"]),
                        "\t".join(["A3", "B3", "C3", "D3", "X3", "Y3", "Z3"]),
                        "\t".join(["A4", "B4", "C4", "D4", "X4", "Y4", "Z4"]),
                        "\t".join(["A5", "B5", "C5", "D5", "X5", "Y5", "Z5"])])
    self.assertEqual(expect + "\n", out_strm.getvalue())

  def test_with_non_matching_header(self):
    # The program should fail if one file has a header and the other doesn't
    # and we don't know the missing value. If the missing value is provided,
    # then it should use this to populate the missing header values in the
    # output.
    pass

  def test_missing_key(self):
    """Program should gracefully handle user giving a key
       name that doesn't exist.
    """
    pass

  def test_duplicate_col_headers(self):
    """If headers are present, each column must have a unique name."""
    pass


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
  _main(sys.argv[1:], sys.argv[0])
