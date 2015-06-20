"""
Date of Creation: 20th June 2015.

Description:    Functions and classes for handling multiple hypothesis testing

Copyright (C) 2010
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

# rpy2 imports
from rpy2.robjects import r

###############################################################################
#                                CONSTANTS                                    #
###############################################################################

ACCEPTED_PVAL_CORRECTION_METHODS = ["BH"]


###############################################################################
#                            P-VALUE CORRECTION                               #
###############################################################################

def correct_pvals(pvals, method="BH", verbose=False):
  """
  Apply multiple hypothesis testing correction to a list of p-values.

  :param pvals: a list or similar iterable containing the p-values to correct;
                can be floats or strings.
  :param method: which correction method to use. Supported options are: BH
  :param verbose: if True, output status messages about the conversion to stderr

  :return: list of corrected p-values
  """
  if verbose:
    sys.stderr.write("Correcting p-values for multiple hypothesis testing...\n")
  pvalstr = ",".join([str(i) for i in pvals])
  r("pvals = c(" + pvalstr + ")")
  return r("p.adjust(pvals, method=\"" + method + "\")")


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################

class TestMHT(unittest.TestCase):

  """Unit tests for the multiple hypothesis testing module."""

  def setUp(self):
    """Set up values for testing the functions in this module."""
    # valid list of p-values, as strings
    self.unif_n_1 = ["0.257441547", "0.606218555", "0.031057321",
                     "0.392350438", "0.089832942", "0.679725728",
                     "0.845886284", "0.904573638", "0.163616912",
                     "0.173740437", "0.882236479", "0.194269205",
                     "0.170241396", "0.690203175", "0.309525257",
                     "0.472407810", "0.190647958", "0.811357538",
                     "0.364202619", "0.974816788", "0.480389155"]
    # valid list of p-values as floats
    self.unif_n_2 = [0.324748977, 0.950240382, 0.427413768, 0.292361543,
                     0.199625459, 0.061821290, 0.487314157, 0.486503664,
                     0.573904823, 0.287293614, 0.883919936, 0.501978056,
                     0.687810603, 0.484497998, 0.012615540, 0.532123089,
                     0.007838910, 0.863231609, 0.313874615, 0.752254812,
                     0.351590089, 0.209198164, 0.706526435, 0.885613515,
                     0.822396869, 0.212189966, 0.444633245, 0.368968403,
                     0.888689927, 0.482703795]

    # corrected self.unif_n_1 by BH method
    self.corrected_unif_n_1_bh = [0.6757841, 0.9058917, 0.5828076, 0.7490327,
                                  0.5828076, 0.9058917, 0.9498023, 0.9498023,
                                  0.5828076, 0.5828076, 0.9498023, 0.5828076,
                                  0.5828076, 0.9058917, 0.7222256, 0.7760133,
                                  0.5828076, 0.9498023, 0.7490327, 0.9748168,
                                  0.7760133]
    # corrected self.unif_n_2 by BH method
    self.corrected_unif_n_2_bh = [0.7925969, 0.9502404, 0.7925969, 0.7925969,
                                  0.7925969, 0.6182129, 0.7925969, 0.7925969,
                                  0.8198640, 0.7925969, 0.9193344, 0.7925969,
                                  0.9193344, 0.7925969, 0.1892331, 0.7981846,
                                  0.1892331, 0.9193344, 0.7925969, 0.9193344,
                                  0.7925969, 0.7925969, 0.9193344, 0.9193344,
                                  0.9193344, 0.7925969, 0.7925969, 0.7925969,
                                  0.9193344, 0.7925969]

  def test_correct_pvals(self):
    """test the correction of p-values for multiple hypothesis testing."""
    output_1 = correct_pvals(self.unif_n_1)
    output_2 = correct_pvals(self.unif_n_2)
    self.assertEqual

    self.assertEqual(len(self.corrected_unif_n_1_bh), len(output_1))
    for i in range(0, len(output_1)):
      self.assertAlmostEqual(self.corrected_unif_n_1_bh[i], output_1[i])
    self.assertEqual(len(self.corrected_unif_n_2_bh), len(output_2))
    for i in range(0, len(output_2)):
      self.assertAlmostEqual(self.corrected_unif_n_2_bh[i], output_2[i])

###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
