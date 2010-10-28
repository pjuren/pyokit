#!/usr/bin/python

""" 
  Date of Creation: 28th May 2010     
                      
  XPIPE Component 
  Description:

  Copyright (C) 2010  
  University of Southern California,
  Philip J. Uren,
  Jin H. Park,
  Andrew D. Smith
  
  Authors: Philip J. Uren, Jin H. Park, Andrew D. Smith
  
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
  
  --------------------
  
  Known Bugs:    None
  
  Revision 
  History:       18th September 2010 -- Philip Uren
                   * added numSequences method 
  
  
  TODO:          None
"""

from collections import deque
from fastaread import FastaRead

def numSequences(fileh):
  """
    DESCRP: Determine how many sequences there are in the file stream passed
            Note: stream will be consumed 
    PARAMS: fileh -- the stream to read the fasta data from, or a string 
                     giving the filename of the file to load 
    RETURN: int -- number of sequences in this fasta stream/file
  """
  if type(fileh).__name__ == "str" : fh = open(fileh)
  else : fh = fileh
  count = 0
  for seq in fastaIterator(fh) : 
    count += 1
  return count

def _isSequenceHead(line):
  """
    DESCRP: Returns true if <line> is the header line of a fastq read sequence
  """
  if len(line) < 1 : return False
  if line[0] == '>' : return True
  return False
      
def fastaIterator(fn, useMutableString = False):
  """
    DESCRP: A generator function which yields reads from <filename>
            the iterator is exhausted when all reads have been 
            read from <filename>
    PARAMS: fn  -  if this is a string, we treat it as a filename, else
                   we treat it as a file-like object, with a readline()
                   method
  """
  prevLine = None
  fh = fn
  if type(fh).__name__ == "str" : fh = open(fh)
  
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
    yield FastaRead(name, seqdata, useMutableString)
    
    # remember where we stopped for next call, or finish 
    prevLine = line
    if prevLine == "" : break