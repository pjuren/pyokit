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
import sys, unittest, os
from collections import deque

# pyokit imports
from pyokit.testing.dummyfiles import DummyInputStream, DummyOutputStream
from pyokit.datastruct.sequence import FastqSequence
from pyokit.util.progressIndicator import ProgressIndicator
from pyokit.util.fileUtils import linesInFile

class FastqFileFormatError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)

def _isSequenceHead(line, prevLine = None):
  """
    @summary Determine whether a string represents a sequence header line in
             fastq format.
    @param line       The line to check
    @param prevLine   @todo
    @return True if <line> is the header line of a fastq read sequence
  """
  # if we know the previous line and it's a quality header, don't call this
  # a sequence header -- the reason is that sanger format uses the '@'
  # symbol in it's quality data and we don't want to mistake that for the
  # a sequence header if it happens to appear as the first character
  #print "check if " + line + " is seq header, with prev line as " + str(prevLine)
  #if prevLine != None and _isQualityHead(prevLine) : return False
  #print "\tpassed prev line is qual header check"

  if len(line) < 1 : return False
  #print "\tpassed empty line test"
  if line[0] == '@' : return True
  #print "\tfirst char not an '@', so it's not a seq header.."
  return False

def _isQualityHead(line):
  """
    @summary Determine whether a string represents a quality header in fastq
             format.
    @return True if <line> is the header line of a fastq quality sequence
  """
  if len(line) < 1 : return False
  if line[0] == '+' : return True
  return False


def fastqIteratorSimple(fn, verbose = False, allowNameMissmatch = False):
  """
    @summary: process a fastq file, but don't take into consideration
              the possibility of data or quality spanning
              more than one line
    @param fn: filename or stream to read data from
    @param ignoreNameMismatch: don't throw error if name in sequence data
                               and quality data parts of a read don't match
    @param verbose: if true, output additional status messages to
                    stderr about progress
  """
  fh = fn
  if type(fh).__name__ == "str" : fh = open(fh)

  # try to get an idea of how much data we have...
  if verbose :
    try :
      totalLines = os.path.getsize(fh.name)
      pind = ProgressIndicator(totalToDo = totalLines,
                               messagePrefix = "completed",
                               messageSuffix = "of processing " +\
                                               fh.name)
    except AttributeError :
      sys.stderr.write("fastqIterator -- warning: " +\
                       "unable to show progress for stream")
      verbose = False

  while True :
    # read four lines.. if we can't get four lines, something is wrong
    lines = []
    gotLines = 0
    while gotLines < 4 :
      l = fh.readline()
      if verbose :
        pind.done = fh.tell()
        pind.showProgress()

      if l == "" :
        # end of file found...
        if gotLines == 0 :
          # ok, not in the middle of a sequence
          break
        else :
          raise FastqFileFormatError("reached end of file in the " +
                                     "middle of sequence data")

      l = l.strip()
      if l == "" : continue
      lines.append(l)
      gotLines += 1

    # couldn't get any more data.. we're done
    if gotLines == 0 : break

    # got our 4 lines, assemble our read..
    # first check that names match
    if lines[0][1:] != lines[2][1:] and not allowNameMissmatch:
      raise FastqFileFormatError("names in sequence don't match : " +\
                                 str(lines[0][1:]) + " != " + str(lines[2][1:]))
    name = lines[0][1:]
    seq = lines[1]
    qual = lines[3]
    yield FastqSequence(name, seq, qual)

def fastqIterator(fn, useMutableString = False, verbose = False, debug = False,
                  sanger = False, allowNameMissmatch = False) :
  """
    @summary
  """
  it = fastqIteratorSimple(fn, verbose=verbose,
                           allowNameMissmatch=allowNameMissmatch)
  for s in it : yield s

def fastqIteratorComplex(fn, useMutableString = False, verbose = False,
                         debug = False, sanger = False):
  """
    @summary: process a fastq file, taking into consideration the
              possibility that sequence or quality data might span
              more than one line
    @param fn: filename or stream to read data from
    @param verbose: if true, output additional status messages to
                    stderr about progress
  """
  fh = fn
  if type(fh).__name__ == "str" : fh = open(fh)
  prevLine = None

  # try to get an idea of how much data we have...
  if verbose :
    try :
      totalLines = linesInFile(fh.name)
      pind = ProgressIndicator(totalToDo = totalLines,
                               messagePrefix = "completed",
                               messageSuffix = "of processing " +\
                                               fh.name)
    except AttributeError :
      sys.stderr.write("fastqIterator -- warning: " +\
                       "unable to show progress for stream")
      verbose = False

  while True :
    # either we have a sequence header left over from the
    # prev call, or we need to read a new one from the file...
    # try to do that now
    seqHeader = ""
    if prevLine != None and not _isSequenceHead(prevLine) :
      raise FastqFileFormatError("terminated on non-read header: " + prevLine)
    if prevLine == None :
      while seqHeader.strip() == "" :
        nl = fh.readline()
        if verbose : pind.done += 1
        if nl == "" :
           ## we're looking for a sequence header, but we're getting eof --
           ## file is empty!
           raise StopIteration()
        seqHeader = nl
    else: seqHeader = prevLine
    name = seqHeader[1:].strip()

    # now we need to read lines until we hit a quality header
    # this is our sequence data
    line = None
    seqdata = ""
    while line == None or not _isQualityHead(line) :
      line = fh.readline()
      if verbose : pind.done += 1
      if line == "" :
        raise FastqFileFormatError("ran out of lines before finding qual head: ")
      if not _isQualityHead(line) : seqdata += line.strip()

    # <line> is now a qual header, keep reading until we see
    # a sequence header.. or we run out of lines.. this is our quality data
    if debug: sys.stderr.write("found quality header: " + line.strip() + "\n")
    prevLine = line
    qualdata = ""
    while not _isSequenceHead(line, prevLine) :
      prevLine = line
      line = fh.readline()
      if verbose : pind.done += 1
      if line == "" : break
      elif not _isSequenceHead(line, prevLine) : qualdata += line.strip()
    if debug:
      sys.stderr.write("finished reading quality data, found: " + qualdata.strip() + "\n")
      sys.stderr.write("loop terminated with line: " + line.strip() + "\n")
      sys.stderr.write("prev when loop termianted was: " + prevLine.strip() + "\n")
    if qualdata.strip() == "" :
      raise FastqFileFormatError("missing quality data..")

    # package it all up..
    yield FastqSequence(name, seqdata, qualdata, useMutableString)
    if verbose : pind.showProgress()

    # remember where we stopped for next call, or finish
    prevLine = line
    if debug : sys.stderr.write("setting prev line to: " + str(prevLine) + "\n")
    if prevLine == "" : break


class FastQUintTests(unittest.TestCase):
  def SetUp(self):
    sanger1 = ()

  def testSangerQual(self):
    """
      @summary: test sequences in sanger format that contain the '@'
                symbol within their quality strings
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
    expect = [FastqSequence(seq1Name, seq1Data, seq1Qual),
              FastqSequence(seq2Name, seq2Data, seq2Qual)]

    ins = DummyInputStream(instr)

    seqs = []
    for seq in fastqIterator(ins):
      seqs.append(seq)

    seqs.sort(key = lambda x: x.sequenceName)
    expect.sort(key = lambda x: x.sequenceName)

    if debug :
      sys.stderr.write("expect\n")
      for seq in expect :
        sys.stderr.write(str(seq) + "\n" + "----------------\n")
      sys.stderr.write("got\n")
      for seq in seqs :
        sys.stderr.write(str(seq) + "\n" + "----------------\n")

    self.assertTrue(seqs == expect)


  def testSangerQual2(self):
    """
      @summary: test sequences in sanger format that contain the '+'
                symbol within their quality strings
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
    expect = [FastqSequence(seq1Name, seq1Data, seq1Qual),
              FastqSequence(seq2Name, seq2Data, seq2Qual)]

    ins = DummyInputStream(instr)

    seqs = []
    for seq in fastqIterator(ins, debug = debug):
      seqs.append(seq)

    seqs.sort(key = lambda x: x.sequenceName)
    expect.sort(key = lambda x: x.sequenceName)

    if debug :
      sys.stderr.write("expect\n")
      for seq in expect :
        sys.stderr.write(str(seq) + "\n" + "----------------\n")
      sys.stderr.write("got\n")
      for seq in seqs :
        sys.stderr.write(str(seq) + "\n" + "----------------\n")

    self.assertTrue(seqs == expect)

if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])
