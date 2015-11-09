"""
Date of Creation: 28th August 2015.

Description:    Join two data tables using a key field in both

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

# standard python imports
import sys
import unittest
import StringIO

# for unit tests
import mock

# pyokit imports
from pyokit.io.fastqIterators import fastqIterator
from pyokit.interface.cli import CLI, Option
from pyokit.util.testing import build_mock_open_side_effect


###############################################################################
#                                  CONSTANTS                                  #
###############################################################################

DEFAULT_VERBOSITY = False


###############################################################################
#                               USER INTERFACE                                #
###############################################################################

def getUI(prog_name, args):
  longDescription = "compute the nucleotide frequency at each position in " +\
                    "a fastq file"
  shortDescription = longDescription

  ui = CLI(prog_name, shortDescription, longDescription)
  ui.minArgs = 1
  ui.maxArgs = 1
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output resultant table to this file. " +
                                  "If omitted, output is to stdout",
                      required=False, type=str))
  ui.addOption(Option(short="h", long="help",
                      description="show this help message ", special=True))
  ui.addOption(Option(short="v", long="verbose",
                      description="outputs additional messages to stdout " +
                                  "about run (default: " +
                                  str(DEFAULT_VERBOSITY) + ")",
                      default=DEFAULT_VERBOSITY, required=False))
  ui.addOption(Option(short="u", long="test",
                      description="run unit tests ", special=True))

  ui.parseCommandLine(args)
  return ui


###############################################################################
#                     COMMAND LINE PROCESSING AND DISPATCH                    #
###############################################################################

def _main(args, prog_name):
  # get options and arguments
  ui = getUI(prog_name, args)

  if ui.optionIsSet("test"):
    # just run unit tests
    unittest.main(argv=[sys.argv[0]])
  elif ui.optionIsSet("help"):
    # just show help
    ui.usage()
  else:
    verbose = (ui.optionIsSet("verbose") is True) or DEFAULT_VERBOSITY

    # get input handles -- we required one args, so we know these will be here.
    infh = open(ui.getArgument(0))

    # make output handle
    outfh = sys.stdout
    if ui.optionIsSet("output"):
      outfh = open(ui.getValue("output"), "w")

    counts = count(infh, verbose)
    write_output(counts, outfh)


###############################################################################
#                             MAIN PROGRAM LOGIC                              #
###############################################################################

def count(fn, verbose=False):
  res = []
  for read in fastqIterator(fn, verbose=verbose, allowNameMissmatch=True):
    for i in range(0, len(read)):
      nuc = read[i]
      if len(res) <= i:
        res.append({})
      if nuc not in res[i]:
        res[i][nuc] = 0
      res[i][nuc] += 1
  return res


def write_output(counts, outfh):
  for i in range(0, len(counts)):
    kys = counts[i].keys()
    kys.sort()
    for nuc in counts[i]:
      outfh.write(str(i) + "\t" + str(nuc) + "\t" + str(counts[i][nuc]) + "\n")


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestNucDist(unittest.TestCase):

  @mock.patch('__builtin__.open')
  def test_simple(self, mock_open):
    out_strm = StringIO.StringIO()
    f = "@700156R:635:C7RW1ACXX:8:1101:1499:2065 2:N:0:41\n" +\
        "ACGAC\n" +\
        "+\n" +\
        "B<<BB\n" +\
        "@700156R:635:C7RW1ACXX:8:1101:1421:2159 2:N:0:41\n" +\
        "TCGAC\n" +\
        "+\n" +\
        "<B0BF\n"

    f_map = {"in.fq": f}
    streams = {"out.dat": out_strm}
    mock_open.side_effect = build_mock_open_side_effect(f_map, streams)

    _main(["-o", "out.dat", "in.fq"], sys.argv[0])
    expect = "\n".join(["\t".join(["0", "A", "1"]),
                        "\t".join(["0", "T", "1"]),
                        "\t".join(["1", "C", "2"]),
                        "\t".join(["2", "G", "2"]),
                        "\t".join(["3", "A", "2"]),
                        "\t".join(["4", "C", "2"])])
    self.assertEqual(expect + "\n", out_strm.getvalue())


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
  try:
    _main(sys.argv[1:], sys.argv[0])
  except Exception as e:
    sys.stderr.write("An unexpected exception occurred. Details: " + str(e) +
                     " \n")
