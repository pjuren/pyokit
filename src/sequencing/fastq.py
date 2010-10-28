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
  History:       18th August 2010 -- Philip Uren
                   * added missing FastqFileFormatError
                   * fixed fastqiterator getting stuck with empty files
                 20th September -- Philip Uren
                   * added mutable string option for fastqIterator 
  
  
  TODO:          None
"""

from collections import deque
from read import FastqRead

class FastqFileFormatError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)

def _isSequenceHead(line):
  """
    DESCRP: Returns true if <line> is the header line of a fastq read sequence
  """
  if len(line) < 1 : return False
  if line[0] == '@' : return True
  return False

def _isQualityHead(line):
  """
    DESCRP: Returns true if <line> is the header line of a fastq quality sequence
  """
  if len(line) < 1 : return False
  if line[0] == '+' : return True
  return False
      
def fastqIterator(fn, useMutableString = False):
  """
    DESCRP: A generator function which yields reads from <filename>
            the iterator is exhausted when all reads have been 
            read from <filename>
    PARAMS: fn  -  if this is a string, we treat it as a filename, else
                   we treat it as a file-like object, with a readline()
                   method
            useMutableString -- if True, underlying sequence objects are
                                constructed as mutable strings, making
                                it possible to rapidly modify them 
                                
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
        nl = fh.readline()
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
      if line == "" : 
        raise FastqFileFormatError("ran out of lines before finding qual head: ")
      if not _isQualityHead(line) : seqdata += line.strip()
      
    # <line> is now a qual header, keep reading until we see 
    # a sequence header.. or we run out of lines.. this is our quality data
    line = None
    qualdata = ""
    while line == None or not _isSequenceHead(line) :
      line = fh.readline()
      if line == "" : break
      elif not _isSequenceHead(line) : qualdata += line.strip()
      
    # package it all up..
    yield FastqRead(name, seqdata, qualdata, useMutableString)
    
    # remember where we stopped for next call, or finish 
    prevLine = line
    if prevLine == "" : break