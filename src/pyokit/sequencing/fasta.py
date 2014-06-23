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
from pyokit.sequencing.fastaread import FastaRead
from pyokit.util.progressIndicator import ProgressIndicator

def numSequences(fileh):
  """
    @summary Determine how many sequences there are in the file stream passed
            Note: stream will be consumed
    @param fileh the stream to read the fasta data from, or a string
                 giving the filename of the file to load. Note that if a stream
                 is given, it will be consumed by this function
    @return the number of sequences in this fasta stream/file
  """
  if type(fileh).__name__ == "str" : fh = open(fileh)
  else : fh = fileh
  count = 0
  for seq in fastaIterator(fh) :
    count += 1
  return count

def _isSequenceHead(line):
  """
    @summary Determine whether a string conforms to the requirements for a
             header line in fasta format.
    @return true if <line> is the header line of a fastq read sequence.
  """
  if len(line) < 1 : return False
  if line[0] == '>' : return True
  return False

def fastaIterator(fn, useMutableString = False, verbose = False):
  """
    @summary A generator function which yields reads from <filename>
             the iterator is exhausted when all reads have been
             read from <filename>
    @param fn       a stream or a string; if this is a string, we treat it as a
                    filename, else we treat it as a file-like object, with a
                    readline() method
    @param verbose  if True, output additional status messages to stderr
                    about progress
  """
  prevLine = None
  fh = fn
  if type(fh).__name__ == "str" : fh = open(fh)

  if verbose :
    try :
      #total = 0
      #for s in fastaIterator(fn.name, verbose = False) : total += 1
      total = os.path.getsize(fh.name)
      pind = ProgressIndicator(totalToDo = total,
                                     messagePrefix = "completed",
                                     messageSuffix = "of processing " +\
                                                      fh.name)
    except AttributeError :
      sys.stderr.write("Warning: " +\
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
        seqHeader = fh.readline()
    else: seqHeader = prevLine
    name = seqHeader[1:].strip()

    # now we need to read lines until we hit another sequence header, or we
    # run out of lines..
    # this is our sequence data

    line = None
    seqdata = ""
    while line == None or not _isSequenceHead(line) :

      line = fh.readline()
      if line == "" :
        # file is finished...
        break
      if not _isSequenceHead(line) : seqdata += line.strip()

    # package it all up..
    if verbose :
      pind.done = fh.tell()
      pind.showProgress()
    yield FastaRead(name, seqdata, useMutableString)

    # remember where we stopped for next call, or finish
    prevLine = line
    if prevLine == "" : break
