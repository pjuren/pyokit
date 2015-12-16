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

# base packes imports
import sys
import os
import unittest

# pyokit imports
from pyokit.testing.dummyfiles import DummyInputStream
from pyokit.datastruct.read import NGSRead
from pyokit.util.progressIndicator import ProgressIndicator
from pyokit.util.fileUtils import linesInFile


class FastqFileFormatError(Exception):
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


def _isSequenceHead(line):
  """
    Determine if a string represents a sequence header line in fastq format.
    Sequence headers start with the '@' symbol, and are followed by sequence
    data.

    :param line:      The line/string to check
    :return: True if the given line conforms to the sequence header format for
             a fastq sequence
  """

  if len(line) < 1:
    return False
  if line[0] == '@':
    return True
  return False


def _isQualityHead(line):
  """
    Determine whether a string represents a quality header in fastq format.
    Quality headers begin with the '+' symbol; subsequent line(s) will be
    quality data.

    :return: True if <line> is the header line of a fastq quality sequence
  """
  if len(line) < 1:
    return False
  if line[0] == '+':
    return True
  return False


def fastqIteratorSimple(fn, verbose=False, allowNameMissmatch=False):
  """
    A generator function that yields FastqSequence objects read from a
    fastq-format stream or filename. This is iterator requires that all
    sequence and quality data is provided on a single line -- put another way,
    it cannot parse fastq files with newline characters interspersed in the
    sequence and/or quality strings. That's probably okay though, as fastq
    files tend not to be formated like that (famous last words..).

    :param fn: filename or stream to read data from.
    :param allowNameMismatch:  don't throw error if name in sequence data
                               and quality data parts of a read don't match.
                               Newer version of CASVA seem to output data like
                               this, probably to save space.
    :param verbose: if True, output additional status messages to stderr about
                    progress.
  """
  fh = fn
  if type(fh).__name__ == "str":
    fh = open(fh)

  # try to get an idea of how much data we have...
  if verbose:
    try:
      totalLines = os.path.getsize(fh.name)
      pind = ProgressIndicator(totalToDo=totalLines,
                               messagePrefix="completed",
                               messageSuffix="of processing "
                                             + fh.name)
    except AttributeError:
      sys.stderr.write("fastqIterator -- warning: " +
                       "unable to show progress for stream")
      verbose = False

  while True:
    # read four lines.. if we can't get four lines, something is wrong
    lines = []
    gotLines = 0
    while gotLines < 4:
      l = fh.readline()
      if verbose:
        pind.done = fh.tell()
        pind.showProgress()

      if l == "":
        # end of file found...
        if gotLines == 0:
          # ok, not in the middle of a sequence
          break
        else:
          raise FastqFileFormatError("reached end of file in the "
                                     + "middle of sequence data")

      l = l.strip()
      if l == "":
        continue
      lines.append(l)
      gotLines += 1

    # couldn't get any more data.. we're done
    if gotLines == 0:
      break

    # got our 4 lines, assemble our read..
    # first check that names match
    if lines[0][1:] != lines[2][1:] and not allowNameMissmatch:
      raise FastqFileFormatError("names in sequence don't match : " +
                                 str(lines[0][1:]) + " != " +
                                 str(lines[2][1:]))
    name = lines[0][1:]
    seq = lines[1]
    qual = lines[3]
    yield NGSRead(seq, name, qual)


def fastqIterator(fn, verbose=False, allowNameMissmatch=False):
  """
    A generator function which yields FastqSequence objects read from a file or
    stream. This is a general function which wraps fastqIteratorSimple. In
    future releases, we may allow dynamic switching of which base iterator is
    used.

    :param fn:                 A file-like stream or a string; if this is a
                               string, it's treated as a filename specifying
                               the location of an input fastq file, else it's
                               treated as a file-like object, which must have a
                               readline() method.
    :param useMustableString:  if True, construct sequences from lists of
                               chars, rather than python string objects, to
                               allow more efficient editing. Use with caution.
    :param verbose:            if True, print messages on progress to stderr.
    :param debug:              if True, print debugging messages to stderr.
    :param sanger:             if True, assume quality scores are in sanger
                               format. Otherwise, assume they're in Illumina
                               format.
    :param allowNameMissmatch: don't throw error if name in sequence data and
                               quality data parts of a read don't match. Newer
                               version of CASVA seem to output data like this,
                               probably to save space.
  """
  it = fastqIteratorSimple(fn, verbose=verbose,
                           allowNameMissmatch=allowNameMissmatch)
  for s in it:
    yield s


def fastq_complex_parse_qual(fh, line, prevLine, verbose=False, pind=None):
  prevLine = line
  qualdata = ""
  while not _isSequenceHead(line):
    prevLine = line
    line = fh.readline()
    if verbose and pind is not None:
      pind.done += 1
    if line == "":
      break
    elif not _isSequenceHead(line):
      qualdata += line.strip()
  if qualdata.strip() == "":
    raise FastqFileFormatError("missing quality data..")
  return qualdata, line, prevLine


def fastq_complex_parse_seq(fh, pind=None, verbose=False):
  line = None
  seqdata = ""
  while line is None or not _isQualityHead(line):
    line = fh.readline()
    if verbose and pind is not None:
      pind.done += 1
    if line == "":
      raise FastqFileFormatError("ran out of lines before finding qual head")
    if not _isQualityHead(line):
      seqdata += line.strip()
  return seqdata, line


def fastq_complex_parse_seq_header(fh, prevLine, pind, verbose):
  seqHeader = ""
  if prevLine is not None and not _isSequenceHead(prevLine):
    raise FastqFileFormatError("terminated on non-read header: " + prevLine)
  if prevLine is None:
    while seqHeader.strip() == "":
      nl = fh.readline()
      if verbose:
        pind.done += 1
      if nl == "":
         # we're looking for a sequence header, but we're getting eof --
         # file is empty!
         raise StopIteration()
      seqHeader = nl
  else:
    seqHeader = prevLine
  return seqHeader[1:].strip(), prevLine


def fastqIteratorComplex(fn, useMutableString=False, verbose=False):
  """
    A generator function which yields FastqSequence objects read from a file or
    stream. This iterator can handle fastq files that have their sequence
    and/or their quality data split across multiple lines (i.e. there are
    newline characters in the sequence and quality strings).

    :param fn:                 A file-like stream or a string; if this is a
                               string, it's treated as a filename specifying
                               the location of an input fastq file, else it's
                               treated as a file-like object, which must have a
                               readline() method.
    :param useMustableString:  if True, construct sequences from lists of
                               chars, rather than python string objects, to
                               allow more efficient editing. Use with caution.
    :param verbose:            if True, print messages on progress to stderr.
    :param debug:              if True, print debugging messages to stderr.
    :param sanger:             if True, assume quality scores are in sanger
                               format. Otherwise, assume they're in Illumina
                               format.
  """
  fh = fn
  if type(fh).__name__ == "str":
    fh = open(fh)
  prevLine = None

  # try to get an idea of how much data we have...
  pind = None
  if verbose:
    try:
      totalLines = linesInFile(fh.name)
      pind = ProgressIndicator(totalToDo=totalLines,
                               messagePrefix="completed",
                               messageSuffix="of processing "
                                             + fh.name)
    except AttributeError:
      sys.stderr.write("fastqIterator -- warning: "
                       + "unable to show progress for stream")
      verbose = False

  while True:
    # either we have a sequence header left over from the
    # prev call, or we need to read a new one from the file...
    # try to do that now
    name, prevLine = fastq_complex_parse_seq_header(fh, prevLine,
                                                    pind, verbose)

    # read lines until we hit a qual header --> this is our sequence data
    seqdata, line = fastq_complex_parse_seq(fh, pind, verbose)

    # <line> is now a qual header, keep reading until we see
    # a sequence header.. or we run out of lines.. this is our quality data
    qualdata, line, prevLine = fastq_complex_parse_qual(fh, line, prevLine,
                                                        verbose, pind)

    # package it all up..
    yield NGSRead(seqdata, name, qualdata, useMutableString)
    if verbose:
      pind.showProgress()

    # remember where we stopped for next call, or finish
    prevLine = line
    if prevLine == "":
      break


class FastQUintTests(unittest.TestCase):
  def testSangerQual(self):
    """
      Unit test for the fastqIterator function -- test that it correctly
      handles sequences in sanger format that contain the '@' symbol within
      their quality strings, including at the start.
    """
    debug = False

    seq1Name = "SRR034466/SRR034466.sra.106 YL_CLIP_2_1_367_452 length=36"
    seq1Data = "GTCATGTGGCCTTCTCAGCAGATTCTTTGTCGTATT"
    seq1Qual = "@HIIIGI:@IIIIIIII.ID;7E;2DA?&71.58@$"
    seq2Name = "SRR034466/SRR034466.sra.107 YL_CLIP_2_1_443_763 length=36"
    seq2Data = "GGATTTCATAGTGATTGTCGTATGCCGTCTTCTTCT"
    seq2Qual = "IIIIIIIIIIIID5IIIIIDIIIIIII)<I%72&4I"

    instr = "@" + seq1Name + "\n" +\
            "" + seq1Data + "\n" +\
            "+" + seq1Name + "\n" +\
            "" + seq1Qual + "\n" +\
            "@" + seq2Name + "\n" +\
            "" + seq2Data + "\n" +\
            "+" + seq2Name + "\n" +\
            "" + seq2Qual + ""
    expect = [NGSRead(seq1Data, seq1Name, seq1Qual),
              NGSRead(seq2Data, seq2Name, seq2Qual)]

    ins = DummyInputStream(instr)

    seqs = []
    for seq in fastqIterator(ins):
      seqs.append(seq)

    seqs.sort(key=lambda x: x.name)
    expect.sort(key=lambda x: x.name)

    if debug:
      sys.stderr.write("expect\n")
      for seq in expect:
        sys.stderr.write(str(seq) + "\n" + "----------------\n")
      sys.stderr.write("got\n")
      for seq in seqs:
        sys.stderr.write(str(seq) + "\n" + "----------------\n")

    self.assertTrue(seqs == expect)

  def testSangerQual2(self):
    """
      Unit test for fastqIterator function -- test that it correctly handles
      sequences in sanger format that contain the '+' symbol within their
      quality strings.
    """
    debug = False

    seq1Name = "SRR034466/SRR034466.sra.1032 YL_CLIP_2_1_865_643 length=36"
    seq1Data = "GACCTCCAAGTAGTCGTAAGCCTACTTCTGCTTGAA"
    seq1Qual = "+II6DEIII/I;-&2>,H%&63$()8205$1D7$(3"
    seq2Name = "SRR034466/SRR034466.sra.1033 YL_CLIP_2_1_845_340 length=36"
    seq2Data = "TAGAGATAGGATTCTGGTGTGTCGTATTCCGTCTTC"
    seq2Qual = "IIIIIIIIIIIII.IID,&@*6;2>HI$91%C(II%"

    instr = "@" + seq1Name + "\n" +\
            "" + seq1Data + "\n" +\
            "+" + seq1Name + "\n" +\
            "" + seq1Qual + "\n" +\
            "@" + seq2Name + "\n" +\
            "" + seq2Data + "\n" +\
            "+" + seq2Name + "\n" +\
            "" + seq2Qual + ""
    expect = [NGSRead(seq1Data, seq1Name, seq1Qual),
              NGSRead(seq2Data, seq2Name, seq2Qual)]

    ins = DummyInputStream(instr)

    seqs = []
    for seq in fastqIterator(ins):
      seqs.append(seq)

    seqs.sort(key=lambda x: x.name)
    expect.sort(key=lambda x: x.name)

    if debug:
      sys.stderr.write("expect\n")
      for seq in expect:
        sys.stderr.write(str(seq) + "\n" + "----------------\n")
      sys.stderr.write("got\n")
      for seq in seqs:
        sys.stderr.write(str(seq) + "\n" + "----------------\n")

    self.assertTrue(seqs == expect)

  def test_complex(self):
    """Test parsing a complex fastq file with line breaks in names and data."""
    debug = False

    seq1Name = "SRR034466/SRR034466.sra.1032 YL_CLIP_2_1_865_643 length=36"
    seq1Data = ("GACCTCCAAGTAGTCGT\n" +
                "AAGCCTACTTCTGCTTGAA")
    seq1Qual = ("III6DEII\n" +
                "I/I;-&2>,H%&\n" +
                "63$()8205$1D7$(3")
    seq2Name = "SRR034466/SRR034466.sra.1033 YL_CLIP_2_1_845_340 length=36"
    seq2Data = "TAGAGATAGGATTCTGGTGTGTCGTATTCCGTCTTC"
    seq2Qual = "IIIIIIIIIIIII.IID,&@*6;2>HI$91%C(II%"

    instr = "\n\n" +\
            "@" + seq1Name + "\n" +\
            "" + seq1Data + "\n" +\
            "+" + seq1Name + "\n" +\
            "" + seq1Qual + "\n" +\
            "@" + seq2Name + "\n" +\
            "" + seq2Data + "\n" +\
            "+" + seq2Name + "\n" +\
            "" + seq2Qual + ""
    expect = [NGSRead(seq1Data.replace('\n', ''), seq1Name,
                      seq1Qual.replace('\n', '')),
              NGSRead(seq2Data, seq2Name, seq2Qual)]

    ins = DummyInputStream(instr)

    seqs = []
    for seq in fastqIteratorComplex(ins):
      seqs.append(seq)

    seqs.sort(key=lambda x: x.name)
    expect.sort(key=lambda x: x.name)

    if debug:
      sys.stderr.write("expect\n")
      for seq in expect:
        sys.stderr.write(str(seq) + "\n" + "----------------\n")
      sys.stderr.write("got\n")
      for seq in seqs:
        sys.stderr.write(str(seq) + "\n" + "----------------\n")

    self.assertTrue(seqs == expect)

if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
