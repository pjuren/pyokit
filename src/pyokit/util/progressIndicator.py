"""
  Date of Creation: 20th May 2010
  Description:      This file defines the ProgressIndicatior class, which is
                    used to report progress on tasks to the user via stderr.

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
