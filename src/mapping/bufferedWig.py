#!/usr/bin/python

""" 
  Date of Creation: 17th November 2010     
                       
  Description: This module defines a class for buffered reading of 
               Wig data. At present, this requires the data to be
               provided in a file, since we might make multiple 
               passes of it.

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
  History:       None
  
  
  TODO:          None  
"""

import sys, os, unittest, math

from mapping.wigIterators import wigIterator, ITERATOR_SORTED_START
from util.fileUtils import linesInFile, openFD
from testing.dummyfiles import DummyInputStream, DummyOutputStream
from datastruct.genomicMap import GenomicMap, GenomicMapError 

DEFAULT_BUFFER_SIZE = 1000000
DEFAULT_CHUNK_SIZE = math.ceil(DEFAULT_BUFFER_SIZE / 100)

class BufferedWig:
  def __init__(self, fn, bufferSize = DEFAULT_BUFFER_SIZE, 
               chunkSize = DEFAULT_CHUNK_SIZE, debug = False):
    """
      @summary: constructor for a BufferedWig object
      @param fn: file descriptor for the underlying wig data. 
                 can't be a stream because multiple passes might be made of
                 this data. Filename is fine. 
      @param bufferSize: the max number of elements to keep in memory at any 
                         given time
      @param chunkSize: the number of elements to flush from the buffer when
                        more space is needed.
      @param debug: output debugging messages if set to True
    """
    self.debug = debug
    if bufferSize != DEFAULT_BUFFER_SIZE and chunkSize == DEFAULT_CHUNK_SIZE:
      chunkSize = math.ceil(bufferSize / float(100))
    
    self.chunkSize = chunkSize
    self.bufferSize = bufferSize
    self.fn = fn
    self.vals = {}
    self.buffer = GenomicMap()
    
  def getItem(self, chrom, loc):
    """
      @summary: get the value of the item on chromosome <chrom> at location
                <loc> -- if it's not in the buffer, it'll be loaded first
      @note: the item is assumed to be only 1 position in length, with and
             end at start + 1
      @param chrom: the chromosome that the desired item is on
      @param loc: the genomic coordinate where the desired item starts
    """
    try :
      return self.buffer.getValue(chrom, loc)
    except GenomicMapError :
      self.loadItem(chrom, loc)
      return self.buffer.getValue(chrom, loc)

  def loadItem(self, chrom, loc):
    """
      @summary: load the item on <chrom> at location <loc> into the buffer;
                Other surrounding items will also be loaded until the buffer
                has been filled, or the end of the stream is met. 
      @param chrom: the chromosome that the desired item is on
      @param loc: the genomic coordinate where the desired item starts
      @note: the item is assumed to be just 1 position long, with an 
             end at start + 1
    """
    # if the buffer is full, make some space
    if len(self.buffer) >= self.bufferSize :
      if self.debug :
        sys.stderr.write("flushing buffer.. removing " +\
                         str(self.chunkSize) + " items \n")
      self.buffer.flush(self.chunkSize)
      if self.debug :
        sys.stderr.write("buffer now looks like: \n" +\
                         str(self.buffer) + "\n")
      
    # find the item we need then keep loading items until we fill up the
    # buffer again or we run out of items
    if self.debug : 
      sys.stderr.write("loading new data... ")
    adding = False
    for item in wigIterator(openFD(self.fn), sortedby = ITERATOR_SORTED_START) :
      if len(self.buffer) >= self.bufferSize : break
      if item.chrom == chrom and item.start == loc : 
        adding = True
      if adding : self.buffer.addValue(item, item.chrom, item.start)
    if self.debug : 
      sys.stderr.write("buffer now looks like: \n" + str(self.buffer) + "\n")