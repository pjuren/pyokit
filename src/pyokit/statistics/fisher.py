#!/usr/bin/python

"""
  Date of Creation: 7th Apr 2012
  Description:      Python wrapper for Fisher's exact test in R

  Copyright (C) 2012-2014
  University of Southern California,
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

import sys, os, unittest, copy

from pyokit.testing.dummyfiles import DummyInputStream, DummyOutputStream
from pyokit.util.progressIndicator import ProgressIndicator
from pyokit.util.fileUtils import linesInFile
from rpy2.robjects import r

def fisherExactTest(a, b, c, d, alternative="two.sided"):
  """
    @summary: ...
    @return: a tuple -- the p-value and odds ratio from Fisher's exact test
  """
  r("m=matrix(c("+str(a)+","+str(b)+","+str(c)+","+str(d)+"), nrow=2)")
  r("res=fisher.test(m, alternative=\"" + alternative +"\")")
  return r("res$p.value")[0], r("res$estimate")[0]


class FisherTests(unittest.TestCase) :
  """
    @summary:
  """

  def testFisherExact(self):
    pval, ratio = fisherExactTest(857,   310487,
                                  43058, 28058896,
                                  alternative="greater")
    self.assertAlmostEqual(pval, 1.2056e-54)
    self.assertAlmostEqual(ratio, 1.79862, 5)


if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])
