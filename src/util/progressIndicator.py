import math, sys

class ProgressIndicator:
  def __init__(self, totalToDo, messagePrefix = None, messageSuffix = None):
    self.total = totalToDo
    self.done = 0
    self.prefix = messagePrefix
    self.suffix = messageSuffix
    self.finished = False
    self.previousMsg = None

    if self.prefix == None : self.prefix = ""
    if self.suffix == None : self.suffix = ""
    
  def showProgress(self):
    if not self.finished :
      percent = math.ceil(100 * self.done / float(self.total))
      
      msg ="\r" + self.prefix + " %d%% " % percent + self.suffix
      if self.previousMsg == None or self.previousMsg != msg :
        sys.stderr.write(msg)
        sys.stderr.flush()
      self.previousMsg = msg    
    
      if percent == 100 : 
        self.finished = True
        sys.stderr.write("\n") 