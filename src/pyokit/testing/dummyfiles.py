#!/usr/bin/python

"""
  Date of Creation: 4th November 2010
  Description:      Provides dummy file streams that behave like real files
                    in most respects, but can be initialised just by providing
                    their content as strings

  Copyright (C) 2010-2014
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

import sys

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

  def reset(self):
    self.current = 0

  def tell(self):
    return self.current

  def seek(self, s):
    self.current = s

  def __iter__(self):
    #return self.lines.__iter__()
    return self

  def next(self):
    self.current += 1
    if self.current > (len(self.lines) - 1) : raise StopIteration
    return self.lines[self.current-1]

  def __eq__(self, o):
    if o == None : return False
    if o == sys.stdin: return True
    return False
  def __ne__(self, o):
    if o == None : return True
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
