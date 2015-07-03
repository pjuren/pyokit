"""
Date of Creation: 2nd July 2014.

Description:   Base exceptions classes for Pyokit.

Copyright (C) 2010-2014
Philip J. Uren,

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


class PyokitError(Exception):

  """The base class for all Pyokit exceptions."""

  def __init__(self, msg):
    """constructor for PyokitExceptions."""
    self.value = msg

  def __str__(self):
    """:return: string representation of this exception."""
    return repr(self.value)
