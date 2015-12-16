#!/usr/bin/python

"""
  Date of Creation: 28th May 2010
  Description:      Functions and iterators for processing fasta files.

  Copyright (C) 2010-2014
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
import sys
import os
import math
import unittest
import StringIO

# pyokit imports
from pyokit.datastruct.sequence import Sequence
from pyokit.util.progressIndicator import ProgressIndicator
from pyokit.util.progressIndicator import ProgressIndicatorError
from pyokit.common.pyokitError import PyokitError
from pyokit.util.testing import build_mock_open_side_effect

# for testing
import mock


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class FastaFileFormatError(PyokitError):
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


###############################################################################
#                              HELPER FUNCTIONS                               #
###############################################################################

def num_sequences(fileh):
  """
    Determine how many sequences there are in fasta file/stream.

    :param fileh: the stream to read the fasta data from, or a string giving
                  the filename of the file to load. Note that if a stream is
                  given, it will be consumed by this function
    :return: the number of sequences in this fasta stream/file
  """
  if type(fileh).__name__ == "str":
    fh = open(fileh)
  else:
    fh = fileh
  count = 0
  for _ in fastaIterator(fh):
    count += 1
  return count


def _isSequenceHead(line):
  """
    Determine whether a string conforms to the requirements for a header line
    in fasta format. Fasta header lines must start with a '>'; empty spaces are
    not stripped, so if your line has any whitespace before '>', it will fail
    this check. There is no check for newline characters either, so strings
    with \n at the end will pass, as will those without it. Multi-line strings
    will also pass.

    :param line: the line to check
    :return: True if the passed line matches the fasta specification for a
             header line, else false.
  """
  if len(line) < 1:
    return False
  if line[0] == '>':
    return True
  return False


###############################################################################
#                             ITERATOR FUNCTIONS                              #
###############################################################################

def __read_seq_header(fh, prev_line):
  """
  Given a file handle/stream, read the next fasta header from that stream.

  :param fh:         filehandle to read from
  :param prev_line:  the previous line read from this filestream
  """
  seq_header = ""
  if prev_line is not None and not _isSequenceHead(prev_line):
    raise FastaFileFormatError("terminated on non-read header: " + prev_line)
  if prev_line is None:
    while seq_header.strip() == "":
      seq_header = fh.readline()
  else:
    seq_header = prev_line
  return seq_header


def __read_seq_data(fh):
  """
  :return: a tuple of the sequence data and the last line that was read from
           the file handle/stream
  """
  line = None
  seq_data = ""
  while line is None or not _isSequenceHead(line):
    line = fh.readline()
    if line == "":
      break  # file is finished...
    if not _isSequenceHead(line):
      seq_data += line.strip()
  return seq_data, line


def __build_progress_indicator(fh):
  try:
    total = os.path.getsize(fh.name)
  except (AttributeError, OSError):
    raise ProgressIndicatorError("Failed to get max size for stream")
  pind = ProgressIndicator(totalToDo=total,
                           messagePrefix="completed",
                           messageSuffix="of processing " + fh.name)
  return pind


def fastaIterator(fn, useMutableString=False, verbose=False):
  """
    A generator function which yields fastaSequence objects from a fasta-format
    file or stream.

    :param fn: a file-like stream or a string; if this is a string, it's
               treated as a filename, else it's treated it as a file-like
               object, which must have a readline() method.
    :param useMustableString: if True, construct sequences from lists of chars,
                              rather than python string objects, to allow
                              more efficient editing. Use with caution.
    :param verbose: if True, output additional status messages to stderr about
                    progress
  """
  fh = fn
  if type(fh).__name__ == "str":
    fh = open(fh)

  if verbose:
    try:
      pind = __build_progress_indicator(fh)
    except ProgressIndicatorError as e:
      sys.stderr.write("Warning: unable to show progress for stream. " +
                       "Reason: " + str(e))
      verbose = False

  prev_line = None
  while True:
    seqHeader = __read_seq_header(fh, prev_line)
    name = seqHeader[1:].strip()
    seq_data, prev_line = __read_seq_data(fh)
    if verbose:
      pind.done = fh.tell()
      pind.showProgress(to_strm=sys.stderr)
    yield Sequence(name, seq_data, useMutableString)

    # remember where we stopped for next call, or finish
    if prev_line == "":
      break


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################

class TestFastaIterators(unittest.TestCase):

  """ Unit tests for fasta iterators. """

  def setUp(self):
    self.file1 = ">s1\n" +\
                 "ACTGATCGATGCGCGATGCTAGTGC\n" +\
                 ">s1\n" +\
                 "ACTGATCGATGCGCGATGCTAGTGC\n" +\
                 ">s1\n" +\
                 "ACTGATCGATGCGCGATGCTAGTGC\n" +\
                 ">s1\n" +\
                 "ACTGATCGATGCGCGATGCTAGTGC\n"
    tmp = StringIO.StringIO(self.file1)
    tmp.seek(0, 2)
    self.file1_size = tmp.tell()

  def test_num_sequences(self):
    fh = StringIO.StringIO(self.file1)
    self.assertEqual(num_sequences(fh), 4)

  def test_fasta_iterator(self):
    """ round-trip for a fasta stream """
    fh = StringIO.StringIO(self.file1)
    res = "\n".join([s.to_fasta_str(include_coords=False)
                     for s in fastaIterator(fh)])
    self.assertEqual(res.strip(), self.file1.strip())

  @mock.patch('os.path.getsize')
  @mock.patch('__builtin__.open')
  def test_verbose(self, mock_open, mock_getsize):
    """ Test progress reporting works properly """
    def mock_get_size_se(*args, **kwargs):
      if (args[0] == "in.fa"):
        return self.file1_size
      else:
        raise ValueError("unknown file: " + args[0])
    outfh = StringIO.StringIO()
    infh = StringIO.StringIO(self.file1)
    f_map = {"in.fa": infh}
    mock_open.side_effect = build_mock_open_side_effect({}, f_map)
    mock_getsize.side_effect = mock_get_size_se
    with mock.patch('sys.stderr', outfh):
      res = ""
      pctngs = []
      for s in fastaIterator("in.fa", verbose=True):
        res += s.to_fasta_str(include_coords=False)
        res += "\n"
        pct = math.ceil(100 * infh.tell() / float(self.file1_size))
        pctngs.append(" %d%% " % pct)
      expect = "\r" + "\r".join(["completed" + m + "of processing in.fa"
                                 for m in pctngs]) + "\n"
      self.assertEqual(outfh.getvalue(), expect)

  @mock.patch('__builtin__.open')
  def test_verbose_stream(self, mock_open):
    """ Test progress reporting fails gracefully on stream of unknown len. """
    outfh = StringIO.StringIO()
    f_map = {"in.fa": self.file1}
    mock_open.side_effect = build_mock_open_side_effect(f_map)
    with mock.patch('sys.stderr', outfh):
      res = "\n".join([s.to_fasta_str(include_coords=False)
                       for s in fastaIterator("in.fa", verbose=True)])
    self.assertEqual("Warning: unable to show progress for stream. Reason:" +
                     " 'Failed to get max size for stream'",
                     outfh.getvalue())
    self.assertEqual(res.strip(), self.file1.strip())


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == '__main__':
    unittest.main()
