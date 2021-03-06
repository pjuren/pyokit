"""
Date of Creation: 1st July 2015.

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
import os
import unittest
import StringIO
import functools

# pyokit imports -- general
from pyokit.interface.cli import CLI, Option
from pyokit.interface.uiError import PyokitUIError
from pyokit.io.indexedFile import IndexedFile
from pyokit.io.indexedFile import IndexError

# pyokit imports -- for indexing repeat-masker alignments
from pyokit.datastruct import multipleAlignment
from pyokit.io.repeatmaskerAlignments import repeat_masker_alignment_iterator

# pyokit imports -- for indexing genome alignments in maf format
from pyokit.io.genomeAlignment import GATestHelper
from pyokit.io.genomeAlignment import genome_alignment_iterator
from pyokit.datastruct.genomeAlignment import JustInTimeGenomeAlignmentBlock
from pyokit.io.genomeAlignment import build_genome_alignment_from_file

# support for enumerations in python 2
from enum import Enum

# for unit tests
import mock


###############################################################################
#                            CONSTANTS AND ENUMS                              #
###############################################################################

DEFAULT_VERBOSITY = False


class FileType(Enum):

  """An enumeration of the possible file types supported by the script."""

  genome_alignment = 1
  reapeat_masker_alignment_file = 2


###############################################################################
#                                 INDEXERS                                    #
###############################################################################

def index_repeatmasker_alignment_by_id(fh, out_fh, vebrose=False):
  """Build an index for a repeat-masker alignment file by repeat-masker ID."""
  def extract_UID(rm_alignment):
    return rm_alignment.meta[multipleAlignment.RM_ID_KEY]

  index = IndexedFile(fh, repeat_masker_alignment_iterator, extract_UID)
  index.write_index(out_fh)


def index_genome_alignment_by_locus(fh, out_fh, verbose=False):
  """Build an index for a genome alig. using coords in ref genome as keys."""
  bound_iter = functools.partial(genome_alignment_iterator,
                                 reference_species="hg19", index_friendly=True)
  hash_func = JustInTimeGenomeAlignmentBlock.build_hash
  idx = IndexedFile(fh, bound_iter, hash_func)
  idx.write_index(out_fh, verbose=verbose)


def get_indexer_by_filetype(fileType):
  """Find the right indexer for a user-specified filetype."""
  try:
    if fileType is FileType.genome_alignment:
      return index_genome_alignment_by_locus
  finally:
    raise ValueError("Unsupported file type: " + str(fileType))


def get_indexer_by_file_extension(ext):
  """Find right indexer using file extension."""
  if ext == ".maf" or ext == "maf":
    return index_genome_alignment_by_locus
  else:
    raise ValueError("Unsupported file type: " + str(ext))


###############################################################################
#                              INDEX LOOKUPS                                  #
###############################################################################

def lookup_genome_alignment_index(index_fh, indexed_fh, out_fh=sys.stdout,
                                  key=None, verbose=False):
  """Load a GA index and its indexed file and extract one or more blocks.

  :param index_fh:   the index file to load. Can be a filename or a
                     stream-like object.
  :param indexed_fh: the file that the index was built for,
  :param key:        A single key, iterable of keys, or None. This key will be
                     used for lookup. If None, user is prompted to enter keys
                     interactively.
  """
  # load the genome alignment as a JIT object
  bound_iter = functools.partial(genome_alignment_iterator,
                                 reference_species="hg19", index_friendly=True)
  hash_func = JustInTimeGenomeAlignmentBlock.build_hash
  idx = IndexedFile(record_iterator=bound_iter, record_hash_function=hash_func)
  idx.read_index(index_fh, indexed_fh)

  if key is None:
    while key is None or key.strip() != "":
      sys.stderr.write("[WAITING FOR KEY ENTRY ON STDIN; " +
                       "END WITH EMPTY LINE]\n")
      key = raw_input()
      # we know keys for genome alignments have tabs as delims, so..
      key = '\t'.join(key.split()).strip()
      if key != "":
        out_fh.write(str(idx[key]) + "\n")
      sys.stderr.write("\n")
  else:
    # we know keys for genome alignments have tabs as delims, so..
    key = '\t'.join(key.split())
    out_fh.write(str(idx[key]) + "\n")


def get_lookup_by_filetype(fileType):
  """Find the right indexer for a user-specified filetype."""
  try:
    if fileType is FileType.genome_alignment:
      return lookup_genome_alignment_index
  finally:
    raise ValueError("Unsupported file type: " + str(fileType))


def get_lookup_by_file_extension(ext):
  """Find right indexer using file extension."""
  if ext == ".maf" or ext == "maf":
    return lookup_genome_alignment_index
  else:
    raise ValueError("Unsupported file type: " + str(ext))


###############################################################################
#                              USER INTERFACE                                 #
###############################################################################

def getUI_build_index(prog_name, args):
  """
  build and return a UI object for the 'build' option.

  :param args: raw arguments to parse
  """
  programName = prog_name
  long_description = "Build an index for one or more files."
  short_description = long_description

  ui = CLI(programName, short_description, long_description)
  ui.minArgs = 0
  ui.maxArgs = -1
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output to given file, else stdout",
                      required=False, type=str))
  ui.addOption(Option(short="t", long="type", argName="filename",
                      description="the type of the file. If missing, " +
                      "the script will try to guess the file type. " +
                      "Supported file types are: " +
                      ", ".join([f.name for f in FileType]),
                      required=False, type=str))
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


def getUI_lookup_index(prog_name, args):
  """
  build and return a UI object for the 'build' option.

  :param args: raw arguments to parse
  """
  programName = prog_name
  long_description = "Build an index for one or more files."
  short_description = long_description

  ui = CLI(programName, short_description, long_description)
  ui.minArgs = 2
  ui.maxArgs = 2
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output to given file, else stdout",
                      required=False, type=str))
  ui.addOption(Option(short="k", long="key", argName="key",
                      description="the key to lookup; if missing, the " +
                      "program runs in interactive mode. ",
                      required=False, type=str))
  ui.addOption(Option(short="t", long="type", argName="filename",
                      description="the type of the file. If missing, " +
                      "the script will try to guess the file type. " +
                      "Supported file types are: " +
                      ", ".join([f.name for f in FileType]),
                      required=False, type=str))
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

def __get_indexer(in_fns, selected_type=None):
  """Determine which indexer to use based on input files and type option."""
  indexer = None
  if selected_type is not None:
    indexer = get_indexer_by_filetype(selected_type)
  else:
    if len(in_fns) == 0:
      raise IndexError("reading from stdin, unable to guess input file " +
                       "type, use -t option to set manually.\n")
    else:
      extension = set([os.path.splitext(f)[1] for f in in_fns])
      assert(len(extension) >= 1)
      if len(extension) > 1:
        raise IndexError("more than one file extension present, unable " +
                         "to get input type, use -t option to set manually.\n")
      else:
        indexer = get_indexer_by_file_extension(list(extension)[0])
  assert(indexer is not None)
  return indexer


def __get_lookup(in_fn, selected_type=None):
  """Determine which lookup func to use based on inpt files and type option."""
  lookup_func = None
  if selected_type is not None:
    lookup_func = get_lookup_by_filetype(selected_type)
  else:
    extension = os.path.splitext(in_fn)[1]
    lookup_func = get_lookup_by_file_extension(extension)
  assert(lookup_func is not None)
  return lookup_func


def isHelpString(s):
  """determine whether an option is a help string."""
  norm = s.strip().lower()
  return norm == "help" or norm == "--help" or norm == "-h"


def main(args=[], prog_name=sys.argv[0]):
  """Figure out which of the two options was used (build, get) and dispatch.

  :param args:      the arguments for this script, as a list of string. Should
                    already have had the script name stripped, so we expect the
                    first item in this list to specify the mode ('build' an
                    index or 'get' elements from an indexed file).
  :param prog_name: the name of the script; taken from command line args by
                    default, but you can change it if you want.
  """
  helpStr = "Build indexes and use them to extract elements from         \n" +\
            "indexed files on-demand.                                    \n" +\
            "                                                            \n" +\
            "------------------------------------------------------------\n" +\
            "                      BUILDING AN INDEX                     \n" +\
            "------------------------------------------------------------\n" +\
            "                                                            \n" +\
            "To build an index for a file called foo.maf, where $ is     \n" +\
            "your prompt:                                                \n" +\
            "                                                            \n" +\
            "$ pyokit index build -o foo.idx foo.maf                     \n" +\
            "                                                            \n" +\
            "for more information about the build command, where $ is    \n" +\
            "your prompt, run:                                           \n" +\
            "                                                            \n" +\
            "$ pyokit index build --help                                 \n" +\
            "                                                            \n" +\
            "------------------------------------------------------------\n" +\
            "                           LOOKUP                           \n" +\
            "------------------------------------------------------------\n" +\
            "To run interactive lookup of values using an index , where $\n" +\
            "is your prompt:                                             \n" +\
            "                                                            \n" +\
            "$ pyokit index lookup foo.maf foo.idx                       \n" +\
            "                                                            \n" +\
            "for more information about the lookup command, where $ is   \n" +\
            "your prompt, run:                                           \n" +\
            "                                                            \n" +\
            "$ pyokit index lookup --help                                \n" +\
            "                                                              "

  if len(args) == 0 or (len(args) == 1 and isHelpString(args[0])):
    sys.stderr.write(helpStr + "\n\n")
  else:
    arg_one = args[0].strip().lower()
    if arg_one == "build":
      main_build_index(args[1:])
    elif arg_one == "lookup":
      main_lookup_index(args[1:])
    else:
      raise PyokitUIError("Unknown option " + arg_one)


def main_lookup_index(args=[], prog_name=sys.argv[0]):
  """
  main entry point for the index script.

  :param args:      the arguments for this script, as a list of string. Should
                    already have had the script name stripped. That is, if
                    there are no args provided, this should be an empty list.
  :param prog_name: the name of the script; taken from command line args by
                    default, but you can change it if you want.
  """
  # get options and arguments
  ui = getUI_lookup_index(prog_name, args)

  # just run unit tests
  if ui.optionIsSet("test"):
    unittest.main(argv=[sys.argv[0]])
    sys.exit()

  # just show help
  if ui.optionIsSet("help"):
    ui.usage()
    sys.exit()
  verbose = (ui.optionIsSet("verbose") is True) or DEFAULT_VERBOSITY

  # required exactly two args, so this should be fine.
  in_indexed, in_index = open(ui.getArgument(0)), open(ui.getArgument(1))

  # Get output handle(s). If only one input file, and no output file, output
  # to standard out.
  out_fh = sys.stdout
  if ui.optionIsSet("output"):
    out_fh = open(ui.getValue("output"), "w")

  # get the lookup key
  key = None
  if ui.optionIsSet("key"):
    key = ui.getValue("key")

  # figure out which indexer to use and then index each file
  op_val = ui.getValue("type") if ui.optionIsSet("type") else None
  lookup_func = __get_lookup(ui.getArgument(0), op_val)

  lookup_func(in_index, in_indexed, out_fh, key, verbose)


def main_build_index(args=[], prog_name=sys.argv[0]):
  """
  main entry point for the index script.

  :param args:      the arguments for this script, as a list of string. Should
                    already have had the script name stripped. That is, if
                    there are no args provided, this should be an empty list.
  :param prog_name: the name of the script; taken from command line args by
                    default, but you can change it if you want.
  """
  # get options and arguments
  ui = getUI_build_index(prog_name, args)

  # just run unit tests
  if ui.optionIsSet("test"):
    unittest.main(argv=[sys.argv[0]])
    sys.exit()

  # just show help
  if ui.optionIsSet("help"):
    ui.usage()
    sys.exit()
  verbose = (ui.optionIsSet("verbose") is True) or DEFAULT_VERBOSITY

  # get input file-handle(s); fall back to stdin if none found.
  in_fhs = [open(x) for x in ui.getAllArguments()]
  if in_fhs == []:
    if sys.stdin.isatty():
      sys.stderr.write("[NO INPUT FILE FOUND; WAITING FOR INPUT FROM STDIN]\n")
    in_fhs = [sys.stdin]

  # Get output handle(s). If only one input file, and no output file, output
  # to standard out.
  out_fhs = [sys.stdout]
  if ui.optionIsSet("output"):
    out_fhs = [open(x, "w") for x in ui.getValue("output").split()]
  if len(in_fhs) != len(out_fhs):
    sys.stderr.write("mismatch between number of input files and output files")
    sys.exit(0)

  # figure out which indexer to use and then index each file
  op_val = ui.getValue("type") if ui.optionIsSet("type") else None
  indexer = __get_indexer(ui.getAllArguments(), op_val)

  for in_fh, out_fh in zip(in_fhs, out_fhs):
    indexer(in_fh, out_fh, verbose)


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################

class TestIndex(unittest.TestCase):

  """Unit tests for the index script."""

  def setUp(self):
    """Set up a few MAF files with whole genome alignments for testing."""
    # these are used for testing indexing of whole genome alignments
    # we'll borrow a genome alignment from the tests in the IO module
    self.b1_hg19 = GATestHelper.b1_hg19
    self.ga_maf1 = GATestHelper.maf1
    self.b2_hg19 = GATestHelper.b2_hg19

  @mock.patch('__builtin__.open')
  def test_index_genome_alig(self, mock_open):
    """Do a full round trip test for a simple whole-genome alignment."""
    in_strm = StringIO.StringIO(self.ga_maf1)
    idx_strm = StringIO.StringIO()

    # replace open with mock
    def open_side_effect(*args, **kwargs):
      if not isinstance(args[0], basestring):
        raise TypeError()
      if args[0] == "one.maf":
        return in_strm
      elif args[0] == "one.idx":
        return idx_strm
      raise IOError("No such file")

    mock_open.side_effect = open_side_effect
    main(["build", "-o", "one.idx", "one.maf"])
    idx_strm.seek(0)
    in_strm.seek(0)
    ga = build_genome_alignment_from_file("one.maf", "hg19", "one.idx")
    b1_res = ga.get_blocks("chr22", 1711, 1720)[0]
    b2_res = ga.get_blocks("chr22", 1770, 1780)[0]
    self.assertEqual(b1_res["hg19.chr22"], self.b1_hg19)
    self.assertEqual(b2_res["hg19.chr22"], self.b2_hg19)

  @mock.patch('__builtin__.open')
  def test_index_lookup(self, mock_open):
    """Test lookup of specific key using full UI."""
    in_strm = StringIO.StringIO(self.ga_maf1)
    idx_strm = StringIO.StringIO()
    out_strm = StringIO.StringIO()

    # replace open with mock
    def open_side_effect(*args, **kwargs):
      if not isinstance(args[0], basestring):
        raise TypeError()
      if args[0] == "one.maf":
        return in_strm
      elif args[0] == "one.idx":
        return idx_strm
      elif args[0] == "out.txt":
        return out_strm
      raise IOError("No such file")

    mock_open.side_effect = open_side_effect

    # build and index in idx_strm
    bound_iter = functools.partial(genome_alignment_iterator,
                                   reference_species="hg19")
    hash_func = JustInTimeGenomeAlignmentBlock.build_hash
    idx = IndexedFile(StringIO.StringIO(self.ga_maf1), bound_iter, hash_func)
    idx.write_index(idx_strm)
    idx_strm.seek(0)

    key = "chr22" + "\t" + "1772" + "\t" + "1825"
    main(["lookup", "-o", "out.txt", "-k", key, "one.maf", "one.idx"])
    self.assertEqual(str(idx[key]).strip(), str(out_strm.getvalue()).strip())


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == "__main__":
    main(sys.argv[1:], sys.argv[0])
