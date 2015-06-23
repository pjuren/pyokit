"""
Date of Creation: 22nd June 2015.

Description:    Classes and functions for computing statistics about streams
                of data as new values arrive.

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
import unittest


###############################################################################
#                       CLASS FOR TRACKING RUNNING MEAN                       #
###############################################################################

class RollingMean(object):

  """Keep track of running mean from a stream of data."""

  def __init__(self):
    """Constructor. See class docstring for parameter details."""
    self._mean = None
    self._vals_added = 0

  def add(self, v):
    """Add a new value."""
    self._vals_added += 1
    if self._mean is None:
      self._mean = v
    self._mean = self._mean + ((v - self._mean) / float(self._vals_added))

  @property
  def mean(self):
    """Get the current mean."""
    return self._mean


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestOnlineStats(unittest.TestCase):

  """Unit tests for this module."""

  def test_rolling_mean(self):
    """Test correct computation of rolling means."""
    r = RollingMean()
    r.add(10)
    self.assertEqual(r.mean, 10)
    r.add(10)
    self.assertEqual(r.mean, 10)
    r.add(5)
    self.assertAlmostEqual(r.mean, 8.3333333)
    r.add(55)
    self.assertAlmostEqual(r.mean, 20)

###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
