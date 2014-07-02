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

import sys, os
from collections import deque
from pyokit.datastruct.sequence import FastaSequence
from pyokit.util.progressIndicator import ProgressIndicator

def numSequences(fileh):
  """
    Determine how many sequences there are in fasta file/stream.

    :param fileh: the stream to read the fasta data from, or a string giving the
                  filename of the file to load. Note that if a stream is given,
                  it will be consumed by this function
    :return: the number of sequences in this fasta stream/file
  """
  if type(fileh).__name__ == "str" : fh = open(fileh)
  else : fh = fileh
  count = 0
  for seq in fastaIterator(fh) :
    count += 1
  return count

def _isSequenceHead(line):
  """
    Determine whether a string conforms to the requirements for a header line in
    fasta format. Fasta header lines must start with a '>'; empty spaces are not
    stripped, so if your line has any whitespace before '>', it will fail this
    check. There is no check for newline characters either, so strings
    with \n at the end will pass, as will those without it. Multi-line strings
    will also pass.

    :param line: the line to check
    :return: True if the passed line matches the fasta specification for a
             header line, else false.
  """
  if len(line) < 1 : return False
  if line[0] == '>' : return True
  return False

def fastaIterator(fn, useMutableString = False, verbose = False):
  """
    A generator function which yields fastaSequence objects from a fasta-format
    file or stream.

    :param fn: a file-like stream or a string; if this is a string, it's treated
               as a filename, else it's treated it as a file-like object, which
               must have a readline() method.
    :param useMustableString: if True, construct sequences from lists of chars,
                              rather than python string objects, to allow
                              more efficient editing. Use with caution.
    :param verbose: if True, output additional status messages to stderr about
                    progress
  """
  prevLine = None
  fh = fn
  if type(fh).__name__ == "str" : fh = open(fh)

  if verbose :
    try :
      total = os.path.getsize(fh.name)
      pind = ProgressIndicator(totalToDo = total,
                               messagePrefix = "completed",
                               messageSuffix = "of processing " +\
                                                fh.name)
    except AttributeError :
      sys.stderr.write("Warning: unable to show progress for stream")
      verbose = False


  while True :
    # either we have a sequence header left over from the prev call, or we need
    # to read a new one from the file... try to do that now
    seqHeader = ""
    if prevLine != None and not _isSequenceHead(prevLine) :
      raise FastqFileFormatError("terminated on non-read header: " + prevLine)
    if prevLine == None :
      while seqHeader.strip() == "" :
        seqHeader = fh.readline()
    else: seqHeader = prevLine
    name = seqHeader[1:].strip()

    # now we need to read lines until we hit another sequence header, or we
    # run out of lines.. this is our sequence data
    line = None
    lineWidth = None
    seqdata = ""
    while line == None or not _isSequenceHead(line) :
      line = fh.readline()
      if line == "" : break  # file is finished...
      if not _isSequenceHead(line) :
        seqdata += line.strip()
        if lineWidth == None :
          lineWidth = len(line.strip())

    # package it all up..
    if verbose :
      pind.done = fh.tell()
      pind.showProgress()
    yield FastaSequence(name, seqdata, lineWidth, useMutableString)

    # remember where we stopped for next call, or finish
    prevLine = line
    if prevLine == "" : break
