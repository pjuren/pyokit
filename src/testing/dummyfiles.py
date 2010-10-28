class DummyInputStream :
  def __init__(self, lines, name = None):
    if type(lines).__name__ == "str" :
      lines = lines.split("\n")
    lines = [x + "\n" for x in lines]
    self.lines = lines
    self.current = 0
    self.length = len(self.lines)
    
    self.name = "none"
    if name != None : self.filename = name
  
  def readline(self):
    if self.current >= self.length : 
      return ""
    self.current += 1 
    return self.lines[self.current-1]
  
  def __iter__(self):
    return self.lines.__iter__()
  
  def __eq__(self, o):
    if o == sys.stdin: return True
    return False
  def __ne__(self, o):
    if o == sys.stdin: return False
    return True
  
  def close(self):
    pass
  
  
  
  
  
  
class DummyOutputStream :
  def __init__(self):
    self.stored = []
    self.prev = ""
    
  def write(self, sth):
    sth = self.prev + sth
    if sth.find("\n") != -1 :
      lines = str(sth).split("\n")
      for line in lines :
        if line == "" : continue 
        self.stored.append(line + "\n")
      self.prev = ""
    else :
      self.prev = sth
      
  def close(self):
    if self.prev != "" :
      self.stored.append(self.prev)
    
  def itemsWritten(self):
    return self.stored
  
  def numItemsWritten(self):
    return len(self.itemsWritten())  
  
  def __str__(self):
    return "\n".join(self.itemsWritten())