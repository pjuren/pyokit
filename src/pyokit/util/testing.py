"""
Date of Creation: November 1st 2015

Description: utilities to assist with testing

Copyright (C) 2015
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

# python imports
import StringIO


def build_mock_open_side_effect(string_d, stream_d):
  """
  Build a mock open side effect using a dictionary of content for the files.

  :param string_d: keys are file names, values are string file contents
  :param stream_d: keys are file names, values are stream of contents
  """
  assert(len(set(string_d.keys()).intersection(set(stream_d.keys()))) == 0)

  def mock_open_side_effect(*args, **kwargs):
    if args[0] in string_d:
      return StringIO.StringIO(string_d[args[0]])
    elif args[0] in stream_d:
      return stream_d[args[0]]
    else:
      raise IOError("No such file: " + args[0] + ". Known these strings " +
                    ", ".join(string_d.keys()) + " and these streams " +
                    ", ".join(stream_d.keys()))
  return mock_open_side_effect
