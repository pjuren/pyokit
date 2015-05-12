#!/usr/bin/python

"""
  Date of Creation: 13th September 2010
  Description:   A class to represent strings that can be modified in-place

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


class MutableString(object) :
  """
    Strings in python are immutable. That brings a number of advantages, but
    one problem is that they are expensive to edit. This class implements a
    string as a list of char, which is cheaper to edit, but cannot be used for
    things like dictionary keys due to it's mutability.
  """
  def __init__(self, strng):
    self.list = [ch for ch in strng]

  def __getitem__(self, indx):
    if isinstance(indx, slice):
      return "".join(self.list.__getitem__(indx))
    return self.list[indx]

  def __setitem__(self, k, v):
    self.list[k] = v

  def __len__(self):
    return len(self.list)

  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return self.__dict__ == other.__dict__
    elif type(other).__name__ == "str":
      return self.__eq__(MutableString(other))
    else:
      return False

  def __ne__(self, other):
      return not self.__eq__(other)

  def __str__(self):
    return "".join(self.list)
