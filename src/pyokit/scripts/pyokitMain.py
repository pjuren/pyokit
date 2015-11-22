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
import unittest

# for testing
import mock

# check for rpy2 -- we need to disable some scripts if it is missing
try:
  from rpy2.robjects import r
  have_functioning_rpy2 = True
except ImportError, e:
  have_functioning_rpy2 = False

# pyokit imports -- exceptions
from pyokit.io.ioError import PyokitIOError
from pyokit.common.pyokitError import PyokitError

# pyokit imports -- scripts
if have_functioning_rpy2:
  from pyokit.scripts import fdr
from pyokit.scripts import index
from pyokit.scripts import conservationProfile
from pyokit.scripts import join
from pyokit.scripts import regionCollapse
from pyokit.scripts import remDupsBED
from pyokit.scripts import convertJunctionReads
from pyokit.scripts import overlapProfile
from pyokit.scripts import readCountExonToGene
from pyokit.scripts import readCounts


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class NoSuchScriptError(PyokitError):
  pass


###############################################################################
#                       DISPATCH HANDLERS FOR SCRIPTS                         #
###############################################################################

def dispatch_fdr(args):
  if have_functioning_rpy2:
    fdr.main(args)
  else:
    sys.stderr.write("The fdr scripts is disabled; it needs rpy2 to " +
                     "run and that package either isn't installed or is " +
                     "not working properly.\n")


def dispatch_index(args):
  index.main(args, "pyokit index")


def dispatch_cons_profile(args):
  conservationProfile.main(args, "pyokit consprofile")


def dispatch_join(args):
  join._main(args, "join")


def dispatch_region_collapse(args):
  regionCollapse.main(args, "pyokit regionCollapse")


def dispatch_convert_junc_reads(args):
  convertJunctionReads._main(args, "remDupsBED")


def dispatch_rem_dups_bed(args):
  remDupsBED._main(args, "remDupsBED")


def dispatch_overlap_profile(args):
  overlapProfile._main(args, "overlapProfile")


def dispatch_read_counts(args):
  readCounts._main(args, "readCounts")


def dispatch_exon_to_gene_read_counts(args):
  readCountExonToGene._main(args, "exonToGeneReadCounts")


dispatchers = {"fdr": dispatch_fdr,
               "index": dispatch_index,
               "consprofile": dispatch_cons_profile,
               "join": dispatch_join,
               "regionCollapse": dispatch_region_collapse,
               "remDupsBED": dispatch_rem_dups_bed,
               "convertJunctionReads": dispatch_convert_junc_reads,
               "overlapProfile": dispatch_overlap_profile,
               "readCounts": dispatch_read_counts,
               "exonToGeneReadCounts": dispatch_exon_to_gene_read_counts}


###############################################################################
#                               MAIN PROG LOGIC                               #
###############################################################################

def dispatch(args):
  """
  Parse the command line and dispatch the appropriate script.

  This function just performs dispatch on the command line that a user
  provided. Basically, look at the first argument, which specifies the
  function the user wants Ribocop to perform, then dispatch to the appropriate
  module.
  """
  prog_name = args[0]
  if prog_name not in dispatchers:
    raise NoSuchScriptError("No such pyokit script: " + prog_name + "\n")
  else:
    dispatchers[prog_name](args[1:])


def main():
  args = sys.argv[1:]
  try:
    dispatch(args)
  except (IOError, PyokitIOError) as e:
    sys.stderr.write("Pyokit - Fatal IOError: " + str(e) + "\n")
    exit(1)
  except PyokitError as e:
    sys.stderr.write("Pyokit - Fatal Error: " + str(e) + "\n")
    exit(1)


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class PyokitSmokeTests(unittest.TestCase):
  """smoke test for pyokit script execution -- just make sure they all run."""

  # @mock.patch('pyokit.scripts.fdr.sys.stderr')
  def test_dispatch(self):
    h_args = ["-h"]
    with mock.patch('sys.stdout'):
      dispatch(["fdr"] + h_args)
      dispatch(["consprofile"] + h_args)
      dispatch(["join"] + h_args)
      dispatch(["regionCollapse"] + h_args)
      dispatch(["remDupsBED"] + h_args)
      dispatch(["convertJunctionReads"] + h_args)
    with mock.patch('sys.stderr'):
      dispatch(["index"] + h_args)


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
