from copy import copy

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
  nfd = copy(fd)
  nfd.reset()
  return nfd

def getFDName(fd):
  if type(fd).__name__ == "str" : return fd
  if type(fd).__name__ == "file" : return fd.name
  return fd.name