#!/usr/bin/python

"""
  Date of Creation: 11th Dec 2014
  Description:      Iterators for processing RepeatMasker files

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
from pyokit.datastruct import retrotransposon


def repeat_masker_iterator(fh, header=True, verbose=False):
  """
  """

  strm = fh
  if type(fh).__name__ == "str":
    strm = open(fh)

  first = True
  for line in strm:
    line = line.strip()
    if line == "":
      continue
    if first and header:
      first = False
      continue
    yield retrotransposon.from_repeat_masker_string(line)
