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

# standard python imports
import sys
import unittest
import inspect

# rpy2 imports -- needed for this module, but we won't demand that the user
# have it -- just disable this module if it cannot be imported
module_disabled = False
try:
  from rpy2.robjects import r
except ImportError:
  sys.stderr.write("Failed to import rpy2; disabling "
                   + inspect.getmodule(inspect.stack()[1][0]).__name__)
  module_disabled = True

if not module_disabled:
  def fisherExactTest(a, b, c, d, alternative="two.sided"):
    """
      @summary: ...
      @return: a tuple -- the p-value and odds ratio from Fisher's exact test
    """
    r("m=matrix(c(" + str(a) + "," + str(b) + "," + str(c) + "," + str(d)
      + "), nrow=2)")
    r("res=fisher.test(m, alternative=\"" + alternative + "\")")
    return r("res$p.value")[0], r("res$estimate")[0]


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class FisherTests(unittest.TestCase) :
  """
    @summary:
  """

  def testFisherExact(self):
    # if the module is disabled, we can't test it... just quietly pass the
    # unit tests
    if module_disabled:
      sys.stderr.write(inspect.getmodule(inspect.stack()[1][0]).__name__
                       + " is disabled, could not import Rpy2; Skipped "
                       + "unit tests.")
    pval, ratio = fisherExactTest(857, 310487,
                                  43058, 28058896,
                                  alternative="greater")
    self.assertAlmostEqual(pval, 1.2056e-54)
    self.assertAlmostEqual(ratio, 1.79862, 5)


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
