
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