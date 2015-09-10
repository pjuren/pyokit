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


class DuplicateKeyError(Exception):

  """..."""

  def __init__(self, msg):
    """..."""
    self.value = msg

  def __str__(self):
    """Get string representation of this exception."""
    return repr(self.value)


class InvalidHeaderError(Exception):

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
    # chomp just the newline char, leave eveerything else alone, so we can
    # handle empty columns in the first and last positions
    line = line.rstrip('\n')

    if verbose:
      pind.done = filehandle.tell()
      pind.showProgress()
    if line == "":
      continue
    yield line


def _build_entry(parts, existing_list_d, key_value, key_field_num,
                 key_is_field_number, header=None, allow_duplicates=False,
                 ignore_missing_keys=False):
  """
  Build and add an entry to existing_list_d.

  If the key is a field number, the entry added will be a list of lists. The
  inner list contains one item per column. and the outer list allows more than
  one entry to be stored per key (but if allow_duplicates is false, an
  exception will be raised if more than one needs to be stored).

  :param parts:               the (already tokenized) list of column entries.
  :param existing_list_d:     a dictionary indexed by key value containing the
                              the already processed entries. The new entry will
                              be added to this using <key_value> as a key
  :param key_value:           the key value to use to add the new entry to
                              <existing_list_d>
  :param key_field_num:       which column is the key_value from (number,
                              indexed from 0)
  :param key_is_field_number: True if the <key_value> is actually the column
                              index, rather than a column name
  :param header:              list giving the names of the columns. Can be None
                              if columns have no names (no header)
  :param allow_duplicates:    ...
  :param ignore_missing_keys: ...
  """
  if key_value.strip() == "":
    if ignore_missing_keys:
      return
    raise MissingKeyError("missing key value")

  if key_value in existing_list_d:
    if not allow_duplicates:
      raise DuplicateKeyError(key_value + " appears multiple times as key")
  else:
    existing_list_d[key_value] = []

  if key_is_field_number:
    # the entry in the dictionary is a list, minus the key field, in the
    # order they occur.
    ne = [parts[i] for i in range(0, len(parts)) if i != key_field_num]
    existing_list_d[key_value].append(ne)
  else:
    # the entry in the dictionary is another dictionary indexed by
    # the header value
    ne = {}
    for i in range(0, len(parts)):
      if i == key_field_num:
        continue
      else:
        ne[header[i]] = parts[i]
    existing_list_d[key_value].append(ne)


def load_file(fn_or_strm, key, key_is_field_number, require_unique_key=True,
              delim="\t", missing_val=None, ignore_missing_keys=False,
              allow_duplicates=False, verbose=False):
  res = {}
  header = None
  key_field_num = key if key_is_field_number else None
  expected_cols_per_line = None

  for line in file_iterator(fn_or_strm, verbose):
    parts = line.split(delim)
    if missing_val is not None:
      parts = [parts[i] if parts[i].strip() != "" else missing_val
               for i in range(0, len(parts))]
    if header is None and not key_is_field_number:
      if any(i.strip() == "" for i in parts):
        raise InvalidHeaderError("header for file has empty fields")
      if len(parts) != len(set(parts)):
        raise InvalidHeaderError("head for file has duplicate names")
      header = parts
      if key not in header:
        raise InvalidHeaderError(key + " does not name a column")
      key_field_num = header.index(key)

      # TODO deal with case where key occurs more than once in the header
    else:
      if expected_cols_per_line is None:
        expected_cols_per_line = len(parts) if header is None else len(header)
      if expected_cols_per_line != len(parts):
        raise ValueError("oops")

      key_val = parts[key_field_num]
      _build_entry(parts, res, key_val, key_field_num, key_is_field_number,
                   header, allow_duplicates=allow_duplicates,
                   ignore_missing_keys=ignore_missing_keys)

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
                                  "Any field that contains only whitespace " +
                                  "is considered missing",
                      required=False, type=str))
  ui.addOption(Option(short="d", long="allow-key-duplicates",
                      description="allow key values to appear more than " +
                                  "once in a file, If this options is set, " +
                                  "all combinations of key-value pairs for " +
                                  "the two files will be output.",
                      required=False))
  ui.addOption(Option(short="i", long="ignore-missing-key",
                      description="skip lines in input files that are " +
                                  "missing a value for the key field ",
                      required=False))
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
  key_one = 0
  key_one_is_field_number = True
  if ui.optionIsSet("field-one"):
    key_one = ui.getValue("field-one")
    try:
      key_one = int(key_one) - 1
      key_one_is_field_number = True
    except ValueError:
      key_one_is_field_number = False

  # get key field in file two
  key_two = 0
  key_two_is_field_number = True
  if ui.optionIsSet("field-two"):
    key_two = ui.getValue("field-two")
    try:
      key_two = int(key_two) - 1
      key_two_is_field_number = True
    except ValueError:
      key_two_is_field_number = False

  # get missing value
  missing_val = None
  if ui.optionIsSet("missing"):
    missing_val = ui.getValue("missing")

  # allow dups?
  allow_duplicate_key_values = ui.optionIsSet("allow-key-duplicates")

  # ignore lines with missing values in the key field?
  ignore_missing_keys = ui.optionIsSet("ignore-missing-key")

  # do our thing..
  process(infh1, infh2, outfh, key_one, key_one_is_field_number, key_two,
          key_two_is_field_number, missing_val, ignore_missing_keys,
          allow_duplicate_key_values, verbose)


###############################################################################
#                             MAIN PROGRAM LOGIC                              #
###############################################################################

def process(infh1, infh2, outfh, key_one, key_one_is_field_number, key_two,
            key_two_is_field_number, missing_val, ignore_missing_keys,
            allow_dup_key_values, verbose=False):
  delim = "\t"
  seen_keys = set()

  mixed_headers = (key_one_is_field_number != key_two_is_field_number)
  if mixed_headers and missing_val is None:
    raise InvalidHeaderError("Cannot join one file with a header to " +
                             "another without one unless you specify " +
                             "a missing-value string for the output " +
                             "(to be placed in the resultant header)")

  f2_dictionary, f2_header = load_file(infh2, key_two, key_two_is_field_number,
                                       missing_val=missing_val,
                                       ignore_missing_keys=ignore_missing_keys,
                                       allow_duplicates=allow_dup_key_values,
                                       verbose=verbose)
  f1_header = None
  key_field_num = key_one if key_one_is_field_number else None
  for line in file_iterator(infh1, verbose):
    parts = line.split(delim)
    if missing_val is not None:
      parts = [parts[i] if parts[i].strip() != "" else missing_val
               for i in range(0, len(parts))]
    if f1_header is None and not key_one_is_field_number:
      if any(i.strip() == "" for i in parts):
        raise InvalidHeaderError("header for first file has empty fields")
      if len(parts) != len(set(parts)):
        raise InvalidHeaderError("header for first file has duplicate names")
      f1_header = parts
      if key_one not in f1_header:
        raise InvalidHeaderError(key_one + " does not name a column")
      key_field_num = f1_header.index(key_one)
      # TODO deal with case where key occurs more than once in the header

      # TODO exception if header is missing and we have no missing val
      # TODO join missing vals otherwise..

      # we know that file one has a header, or we wouldn't be here...
      if not key_two_is_field_number:
        outfh.write(delim.join(f1_header) + delim +
                    delim.join(f2_header) + "\n")
      else:
        dummy_h = len(f2_dictionary[f2_dictionary.keys()[0]][0]) *\
            [missing_val]
        outfh.write(delim.join(f1_header) + delim +
                    delim.join(dummy_h) + "\n")
    else:
      key_val = parts[key_field_num]
      if key_val.strip() == "":
        if ignore_missing_keys:
          continue
        raise MissingKeyError("missing key value")
      if key_val in seen_keys and not allow_dup_key_values:
        raise DuplicateKeyError(key_val + " appears multiple times as key")
      seen_keys.add(key_val)

      if f2_header is not None:
        # f2_dictionary entry is another dictionary indexed by header value
        if key_val not in f2_dictionary:
          continue
        for ne in f2_dictionary[key_val]:
          outfh.write(delim.join(parts) + delim)
          first = True
          for h_val in f2_header:
            if h_val == key_two:
              continue

            if first:
              first = False
            else:
              outfh.write(delim)
            outfh.write(ne[h_val])
          outfh.write("\n")
      else:
        # f1_dictionary entry is a list
        for ne in f2_dictionary[key_val]:
          outfh.write(delim.join(parts) + delim + delim.join(ne) + "\n")


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

def build_mock_open_side_effect(string_d, stream_d):
  """
  Build a mock open side effect using a dictionary of content for the files.

  :param string_d: keys are file names, values are string file contents
  :param stream_d: keys are file names, values are stream of contents
  """
  assert(len(set(string_d.keys()).intersection(set(stream_d.keys()))) == 0)

  def mock_open_side_effect(*args, **kwargs):
    if args[0] in string_d:
      return StringIO.StringIO(string_d[args[0]])
    elif args[0] in stream_d:
      return stream_d[args[0]]
    else:
      raise IOError("No such file: " + args[0])
  return mock_open_side_effect


class TestJoin(unittest.TestCase):

  """Unit tests for this script."""

  def setUp(self):
    # matches cleanly with test_file_two
    self.test_headr_one = "\t".join(["AA", "BB", "CC", "DD"])
    self.test_headr_one_dup = "\t".join(["AA", "BB", "AA", "DD"])
    self.test_file_one = ["\t".join(["A1", "B1", "C1", "D1"]),
                          "\t".join(["A2", "B2", "C2", "D2"]),
                          "\t".join(["A3", "B3", "C3", "D3"]),
                          "\t".join(["A4", "B4", "C4", "D4"]),
                          "\t".join(["A5", "B5", "C5", "D5"])]

    # same as test_file_one, but with empty fields
    self.test_headr_one_gapped = "\t".join(["AA", "BB", "", "DD"])
    self.test_file_one_gapped = ["\t".join(["A1", "B1", "C1", "D1"]),
                                 "\t".join(["A2", "", "C2", "D2"]),
                                 "\t".join(["A3", "", "C3", "D3"]),
                                 "\t".join(["", "B4", "C4", "D4"]),
                                 "\t".join(["  ", "B5", "C5", "D5"])]

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

    self.strs = {"one.dat": "\n".join(self.test_file_one),
                 "two.dat": "\n".join(self.test_file_two),
                 "one_hdr.dat": (self.test_headr_one + "\n" +
                                 "\n".join(self.test_file_one)),
                 "two_hdr.dat": (self.test_headr_two + "\n" +
                                 "\n".join(self.test_file_two)),
                 "one_gpd.dat": "\n".join(self.test_file_one_gapped),
                 "one_hdr_gpd.dat": (self.test_headr_one_gapped + "\n" +
                                     "\n".join(self.test_file_one_gapped)),
                 "one_dup.dat": (self.test_headr_one_dup + "\n" +
                                 "\n".join(self.test_file_one)),
                 "one_dup_fields.dat": (self.test_headr_four + "\n" +
                                        "\n".join(self.test_file_four))}

  @mock.patch('__builtin__.open')
  def test_simple_headerless_join(self, mock_open):
    """Test joining two files where all keys are present in both, no header."""
    out_strm = StringIO.StringIO()
    strings = {"one.dat": "\n".join(self.test_file_one),
               "two.dat": "\n".join(self.test_file_two)}
    streams = {"out.dat": out_strm}
    mock_open.side_effect = build_mock_open_side_effect(strings, streams)

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
    streams = {"out.dat": out_strm}
    mock_open.side_effect = build_mock_open_side_effect(self.strs, streams)

    _main(["-a", "BB", "-b", "BX", "-o", "out.dat",
           "one_hdr.dat", "two_hdr.dat"], "join")
    expect = "\t".join(["AA", "BB", "CC", "DD", "XX", "YY", "ZZ"]) + "\n" +\
             "\n".join(["\t".join(["A1", "B1", "C1", "D1", "X1", "Y1", "Z1"]),
                        "\t".join(["A2", "B2", "C2", "D2", "X2", "Y2", "Z2"]),
                        "\t".join(["A3", "B3", "C3", "D3", "X3", "Y3", "Z3"]),
                        "\t".join(["A4", "B4", "C4", "D4", "X4", "Y4", "Z4"]),
                        "\t".join(["A5", "B5", "C5", "D5", "X5", "Y5", "Z5"])])
    self.assertEqual(expect + "\n", out_strm.getvalue())

  @mock.patch('__builtin__.open')
  def test_join_header_non_header(self, mock_open):
    """Test joining a file with a header to a file without one."""
    out_strm = StringIO.StringIO()
    streams = {"out.dat": out_strm}
    mock_open.side_effect = build_mock_open_side_effect(self.strs, streams)

    # fails if the missing value is not provided...
    args = ["-a", "BB", "-b", "2", "-o", "out.dat", "one_hdr.dat", "two.dat"]
    self.assertRaises(InvalidHeaderError, _main, args, "join")

    # okay if it is.
    _main(["-a", "BB", "-b", "2", "-o", "out.dat", "-m",
           "UNKNOWN", "one_hdr.dat", "two.dat"],
          "join")
    expect = "\t".join(["AA", "BB", "CC", "DD",
                        "UNKNOWN", "UNKNOWN", "UNKNOWN"]) + "\n" +\
             "\n".join(["\t".join(["A1", "B1", "C1", "D1", "X1", "Y1", "Z1"]),
                        "\t".join(["A2", "B2", "C2", "D2", "X2", "Y2", "Z2"]),
                        "\t".join(["A3", "B3", "C3", "D3", "X3", "Y3", "Z3"]),
                        "\t".join(["A4", "B4", "C4", "D4", "X4", "Y4", "Z4"]),
                        "\t".join(["A5", "B5", "C5", "D5", "X5", "Y5", "Z5"])])
    self.assertEqual(expect + "\n", out_strm.getvalue())

  @mock.patch('__builtin__.open')
  def test_default_keys(self, mock_open):
    """Test using the default keys (first field in each file)."""
    out_strm = StringIO.StringIO()
    streams = {"out.dat": out_strm}
    mock_open.side_effect = build_mock_open_side_effect(self.strs, streams)

    _main(["-o", "out.dat", "one.dat", "one.dat"], "join")
    expect = ["\t".join(["A1", "B1", "C1", "D1", "B1", "C1", "D1"]),
              "\t".join(["A2", "B2", "C2", "D2", "B2", "C2", "D2"]),
              "\t".join(["A3", "B3", "C3", "D3", "B3", "C3", "D3"]),
              "\t".join(["A4", "B4", "C4", "D4", "B4", "C4", "D4"]),
              "\t".join(["A5", "B5", "C5", "D5", "B5", "C5", "D5"])]
    self.assertEqual("\n".join(expect) + "\n", out_strm.getvalue())

  @mock.patch('__builtin__.open')
  def test_empty_fields(self, mock_open):
    """Test joining files where there are empty fields."""
    out_strm = StringIO.StringIO()
    streams = {"out.dat": out_strm}
    mock_open.side_effect = build_mock_open_side_effect(self.strs, streams)

    # should fail gracefully when there is an empty field in the header
    # of either file, or both
    args = ["-o", "out.dat", "-a", "BB", "one_hdr.dat", "one_hdr_gpd.dat"]
    self.assertRaises(InvalidHeaderError, _main, args, "join")
    args = ["-o", "out.dat", "-b", "BB", "one_hdr_gpd.dat", "one_hdr.dat"]
    self.assertRaises(InvalidHeaderError, _main, args, "join")
    args = ["-o", "out.dat", "-a", "BB", "-b", "BB",
            "one_hdr_gpd.dat", "one_hdr_gpd.dat"]
    self.assertRaises(InvalidHeaderError, _main, args, "join")

    # should fail gracefully when there is an empty field in the key column
    # unless the option to ignore those rows is given.
    args = ["-o", "out.dat", "-a", "2", "-b", "2", "one_gpd.dat", "one.dat"]
    self.assertRaises(MissingKeyError, _main, args, "join")
    args = ["-o", "out.dat", "-a", "2", "-b", "2", "one.dat", "one_gpd.dat"]
    self.assertRaises(MissingKeyError, _main, args, "join")
    args = ["-o", "out.dat", "-a", "2", "-b", "2", "one_gpd.dat",
            "one_gpd.dat"]
    self.assertRaises(MissingKeyError, _main, args, "join")
    out_strm.truncate(0)
    out_strm.seek(0)
    _main(["-i", "-o", "out.dat", "-a", "2", "-b", "2", "one_gpd.dat",
           "one_gpd.dat"], "join")
    expect = ["\t".join(["A1", "B1", "C1", "D1", "A1", "C1", "D1"]),
              "\t".join(["", "B4", "C4", "D4", "", "C4", "D4"]),
              "\t".join(["  ", "B5", "C5", "D5", "  ", "C5", "D5"])]
    self.assertEqual("\n".join(expect) + "\n", out_strm.getvalue())

    # should work when gaps are not in the key field, regardless of whether
    # missing value is provided or not
    out_strm.truncate(0)
    out_strm.seek(0)
    _main(["-o", "out.dat", "-a", "3", "-b", "3", "one_gpd.dat",
           "one_gpd.dat"], "join")
    expect = ["\t".join(["A1", "B1", "C1", "D1", "A1", "B1", "D1"]),
              "\t".join(["A2", "", "C2", "D2", "A2", "", "D2"]),
              "\t".join(["A3", "", "C3", "D3", "A3", "", "D3"]),
              "\t".join(["", "B4", "C4", "D4", "", "B4", "D4"]),
              "\t".join(["  ", "B5", "C5", "D5", "  ", "B5", "D5"])]
    self.assertEqual("\n".join(expect) + "\n", out_strm.getvalue())
    out_strm.truncate(0)
    out_strm.seek(0)
    _main(["-o", "out.dat", "-a", "3", "-b", "3", "-m", "UNKNOWN",
           "one_gpd.dat", "one_gpd.dat"], "join")
    expect = ["\t".join(["A1", "B1", "C1", "D1", "A1", "B1", "D1"]),
              "\t".join(["A2", "UNKNOWN", "C2", "D2", "A2", "UNKNOWN", "D2"]),
              "\t".join(["A3", "UNKNOWN", "C3", "D3", "A3", "UNKNOWN", "D3"]),
              "\t".join(["UNKNOWN", "B4", "C4", "D4", "UNKNOWN", "B4", "D4"]),
              "\t".join(["UNKNOWN", "B5", "C5", "D5", "UNKNOWN", "B5", "D5"])]
    self.assertEqual("\n".join(expect) + "\n", out_strm.getvalue())

  @mock.patch('__builtin__.open')
  def test_missing_key(self, mock_open):
    """Should gracefully handle user giving a key name that doesn't exist."""
    out_strm = StringIO.StringIO()
    streams = {"out.dat": out_strm}
    mock_open.side_effect = build_mock_open_side_effect(self.strs, streams)

    # should be fine....
    try:
      _main(["-o", "out.dat", "-a", "BB", "-b", "BX", "one_hdr.dat",
             "two_hdr.dat"], "join")
    except Exception:
      self.fail()
    # should fail ..
    args = ["-o", "out.dat", "-a", "BX", "-b", "BX", "one_hdr.dat",
            "two_hdr.dat"]
    self.assertRaises(InvalidHeaderError, _main, args, "join")
    args = ["-o", "out.dat", "-a", "BB", "-b", "BB", "one_hdr.dat",
            "two_hdr.dat"]
    self.assertRaises(InvalidHeaderError, _main, args, "join")

  @mock.patch('__builtin__.open')
  def test_duplicate_col_headers(self, mock_open):
    """If headers are present, each column must have a unique name."""
    out_strm = StringIO.StringIO()
    streams = {"out.dat": out_strm}
    mock_open.side_effect = build_mock_open_side_effect(self.strs, streams)
    args = ["-o", "out.dat", "-a", "BB", "-b", "BX", "one_dup.dat",
            "two_hdr.dat"]
    self.assertRaises(InvalidHeaderError, _main, args, "join")
    args = ["-o", "out.dat", "-a", "BX", "-b", "BB", "two_hdr.dat",
            "one_dup.dat"]
    self.assertRaises(InvalidHeaderError, _main, args, "join")

  def test_failure_on_ragged_data_frame(self):
    # number of elements should be the same on each line...
    # TODO
    pass

  @mock.patch('__builtin__.open')
  def test_fail_on_dup_key_vals(self, mock_open):
    """Duplicate key vals should cause failure when approp. option not set."""
    out_strm = StringIO.StringIO()
    streams = {"out.dat": out_strm}
    mock_open.side_effect = build_mock_open_side_effect(self.strs, streams)
    args = ["-o", "out.dat", "-a", "BB", "-b", "BX", "one_dup_fields.dat",
            "two_hdr.dat"]
    self.assertRaises(DuplicateKeyError, _main, args, "join")
    args = ["-o", "out.dat", "-a", "BX", "-b", "BB", "two_hdr.dat",
            "one_dup_fields.dat"]
    self.assertRaises(DuplicateKeyError, _main, args, "join")

  @mock.patch('__builtin__.open')
  def test_duplicate_key_value_all_combinations(self, mock_open):
    # if a key value appears more than once, by default the program will just
    # exit with an error, but there is an option to process all combinations
    out_strm = StringIO.StringIO()
    streams = {"out.dat": out_strm}
    mock_open.side_effect = build_mock_open_side_effect(self.strs, streams)

    _main(["-d", "-o", "out.dat", "-a", "BB", "-b", "BX", "one_dup_fields.dat",
           "two_hdr.dat"], "join")
    expect = ["\t".join(["AA", "BB", "CC", "DD", "XX", "YY", "ZZ"]),
              "\t".join(["A1", "B1", "C1", "D1", "X1", "Y1", "Z1"]),
              "\t".join(["A2", "B1", "C2", "D2", "X1", "Y1", "Z1"]),
              "\t".join(["A3", "B3", "C3", "D3", "X3", "Y3", "Z3"]),
              "\t".join(["A4", "B3", "C4", "D4", "X3", "Y3", "Z3"]),
              "\t".join(["A5", "B5", "C5", "D5", "X5", "Y5", "Z5"])]
    print "\n".join(expect) + "\n"
    print "**********"
    print out_strm.getvalue()
    self.assertEqual("\n".join(expect) + "\n", out_strm.getvalue())

    out_strm.truncate(0)
    out_strm.seek(0)
    _main(["-d", "-o", "out.dat", "-a", "BX", "-b", "BB", "two_hdr.dat",
           "one_dup_fields.dat"], "join")
    expect = ["\t".join(["XX", "BX", "YY", "ZZ", "AA", "CC", "DD"]),
              "\t".join(["X1", "B1", "Y1", "Z1", "A1", "C1", "D1"]),
              "\t".join(["X1", "B1", "Y1", "Z1", "A2", "C2", "D2"]),
              "\t".join(["X3", "B3", "Y3", "Z3", "A3", "C3", "D3"]),
              "\t".join(["X3", "B3", "Y3", "Z3", "A4", "C4", "D4"]),
              "\t".join(["X5", "B5", "Y5", "Z5", "A5", "C5", "D5"])]
    print "\n".join(expect) + "\n"
    print "**********"
    print out_strm.getvalue()
    self.assertEqual("\n".join(expect) + "\n", out_strm.getvalue())

  def test_duplicate_key_value_join_fields(self):
    # if a key value appears more than once, by default the program will just
    # exit with an error, but there is an option to join the mismatched fields
    # TODO
    pass

  def test_output_unmatched_keys(self):
    # there should be an option to output lines with key values that don't
    # match anything in the other file. This should require that the missing
    # value is provided, to fill in unknown columns
    # TODO
    pass


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
  _main(sys.argv[1:], sys.argv[0])