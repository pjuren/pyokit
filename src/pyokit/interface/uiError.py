#!/usr/bin/python

"""
Date of Creation: 3rd July 2015.

Description:   Generic user interface error exception classes

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

# pyokit imports
from pyokit.common.pyokitError import PyokitError


###############################################################################
#                            IO EXCEPTION CLASSES                             #
###############################################################################

class PyokitUIError(PyokitError):

  """Generic pyokit UI error class."""

  def __init__(self, msg):
    """Constructor for UI errors."""
    self.value = msg

  def __str__(self):
    """:return: string representation of this UI error."""
    return repr(self.value)
