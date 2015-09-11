"""
Date of Creation: 18th June 2015.

Description:   main dispatch for pyokit scripts

Copyright (C) 2010-2015
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

# standard python imports
import sys

# pyokit imports -- exceptions
from pyokit.io.ioError import PyokitIOError
from pyokit.common.pyokitError import PyokitError

# pyokit imports -- scripts
from pyokit.scripts import fdr
from pyokit.scripts import index
from pyokit.scripts import conservationProfile
from pyokit.scripts import join
from pyokit.scripts import regionCollapse


def main():
  """
  Parse the command line and dispatch the appropriate script.

  This function just performs dispatch on the command line that a user
  provided. Basically, look at the first argument, which specifies the
  function the user wants Ribocop to perform, then dispatch to the appropriate
  module.
  """
  try:
    if sys.argv[1] == "fdr":
      fdr.main(sys.argv[2:])
    elif sys.argv[1] == "index":
      index.main(sys.argv[2:], "pyokit index")
    elif sys.argv[1] == "consprofile":
      conservationProfile.main(sys.argv[2:], "pyokit consprofile")
    elif sys.argv[1] == "join":
      join._main(sys.argv[2:], "join")
    elif sys.argv[1] == "regionCollapse":
      regionCollapse.main(sys.argv[2:], "pyokit regionCollapse")
    else:
      sys.stderr.write("Pyokit: I don't recognise the option '" + sys.argv[1] +
                       "'.\n")
  except (IOError, PyokitIOError) as e:
    sys.stderr.write("Pyokit - Fatal IOError: " + str(e) + "\n")
    exit(1)
  except PyokitError as e:
    sys.stderr.write("Pyokit - Fatal Error: " + str(e) + "\n")
    exit(1)
