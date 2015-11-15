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

# standard python imports
import math
import sys
import unittest
import StringIO

# pyokit imports
from pyokit.common.pyokitError import PyokitError


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################
class ProgressIndicatorError(PyokitError):
  pass


###############################################################################
#                           ProgressIndicator CLASS                           #
###############################################################################

class ProgressIndicator(object):
  def __init__(self, totalToDo, messagePrefix=None, messageSuffix=None):
    self.total = totalToDo
    self.done = 0
    self.prefix = messagePrefix
    self.suffix = messageSuffix
    self.finished = False
    self.previousMsg = None

    if self.prefix is None:
      self.prefix = ""
    if self.suffix is None:
      self.suffix = ""

  def showProgress(self, to_strm=sys.stderr):
    if not self.finished:
      percent = math.ceil(100 * self.done / float(self.total))

      msg = "\r" + self.prefix + " %d%% " % percent + self.suffix
      if self.previousMsg is None or self.previousMsg != msg:
        to_strm.write(msg)
        to_strm.flush()
      self.previousMsg = msg

      if percent == 100:
        self.finished = True
        to_strm.write("\n")


###############################################################################
#                                 UNIT TESTS                                  #
###############################################################################

class TestProgressIndicator(unittest.TestCase):

  def test_basic(self):
    """Test that no more than 100 msgs are ever output, and they're right."""
    MAX_VAL = 2000
    outfh = StringIO.StringIO()
    nums = range(0, MAX_VAL)
    pind = ProgressIndicator(totalToDo=len(nums), messagePrefix="done",
                             messageSuffix="of numbers")
    for x in nums:
      pind.done += 1
      pind.showProgress(outfh)

    expect = "\r" + "\r".join(["done " + str(i) + "% of numbers"
                               for i in range(1, 101)])
    expect += "\n"
    self.assertEqual(expect, outfh.getvalue())


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
