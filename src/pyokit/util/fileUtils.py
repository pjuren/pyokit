"""
  Date of Creation: 22nd May 2011
  Description:      functions for manupulating files, filehandles and filenames 

  Copyright (C) 2011-2014
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


import sys, os
from copy import copy
from pyokit.util.progressIndicator import ProgressIndicator

def genericFileIterator(fn, verbose=False):
  """
    @summary: iterate over a file, returning non-blank lines
    @param fn: either a string representing the name of the file or a file object
    @param verbose: if True, output status messages to stderr
  """
  if type(fn).__name__ == "str" : fh = open(fn)
  else : fh = fn

  if verbose :
    try :
      pind = ProgressIndicator(totalToDo = os.path.getsize(fh.name),
                               messagePrefix = "completed",
                               messageSuffix = "of processing " +\
                                               fh.name)
      junk = fh.tell()
    except :
      sys.stderr.write("Cannot show progress for stream.. doesn't behave like a file")
      verbose = False

  for line in fh :
    if verbose :
      pind.done = fh.tell()
      pind.showProgress()
    line = line.strip()
    if line == "" : continue
    yield line

def getUniqueFilename(dir = None, base = None):
  """
    DESCRP: Generate a filename in the directory <dir> which is
            unique (i.e. not in use at the moment)
    PARAMS: dir  -- the directory to look in. If None, use CWD
            base -- use this as the base name for the filename
    RETURN: string -- the filename generated
  """
  while True :
    fn = str(random.randint(0,100000)) + ".tmp"
    if not os.path.exists(fn) : break
  return fn

def linesInFile(fd):
  if type(fd).__name__ == "str" : f = open(fd)
  else : f = fd
  t = sum(1 for line in f)
  f.close()
  return t

def openFD(fd):
  """
    @summary: given a descriptor for a file (e.g. path or file stream), will
              attempt to return a new stream for the file ready to be read
              from the beginning
    @param fd: the file descriptor to open a stream to
    @raise IOError: if a new stream reset to the begining of the file cannot
                    be constructed
  """
  if type(fd).__name__ == "str" : return open(fd)
  if type(fd).__name__ == "file" : return open(fd.name)
  if type(fd).__name__ == "mmap" : return fd
  nfd = copy(fd)
  nfd.reset()
  return nfd

def getFDName(fd):
  if type(fd).__name__ == "str" : return fd
  if type(fd).__name__ == "file" : return fd.name
  return fd.name
