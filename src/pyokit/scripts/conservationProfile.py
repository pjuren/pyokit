"""
Date of Creation: 21st June 2015.

Description:    Given a set of genomic intervals in some genome in BED format
                and a genome alignment with that species as reference, produce
                a profile of conservation centered on those locations.

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
import os
import sys
import copy
import unittest

# pyokit imports
from pyokit.interface.cli import CLI
from pyokit.interface.cli import Option
from pyokit.io.bedIterators import BEDIterator
from pyokit.datastruct.genomeAlignment import NoSuchAlignmentColumnError
from pyokit.datastruct.genomeAlignment import NoUniqueColumnError


###############################################################################
#                                 CONSTANTS                                   #
###############################################################################

# the following are used for specifying which part of BED intervals to use
# as the center of the profile
START = "START"
END = "END"
FIVE_PRIME = "5PRIME"
THREE_PRIME = "3PRIME"
CENTRE = "CENTRE"
WINDOW_CENTRE_OPTIONS = [START, END, CENTRE, THREE_PRIME, FIVE_PRIME]

# default parameter values
DEFAULT_WINDOW_SIZE = 30
DEFAULT_VERBOSITY = False
DEFAULT_WINDOW_CENTRE = CENTRE


###############################################################################
#              HELPER FUNCTIONS FOR TRANSFORMING THE INPUT REGIONS            #
###############################################################################

def center_start(r, window_size):
  """
  Center a region on its start and expand it to window_size bases.

  :return: the new region.
  """
  res = copy.copy(r)
  res.end = res.start + window_size / 2
  res.start = res.end - window_size
  return res


def center_end(r, window_size):
  """
  Center a region on its end and expand it to window_size bases

  :return: the new region.
  """
  res = copy.copy(r)
  res.start = res.end - window_size / 2
  res.end = res.start + window_size
  return res


def center_middle(r, window_size):
  """
  Center a region on its middle and expand it to window_size bases

  :return: the new region.
  """
  res = copy.copy(r)
  mid = res.start + (len(res) / 2)
  res.start = mid - (window_size / 2)
  res.end = res.start + window_size
  return res


def transform_locus(region, window_center, window_size):
  """
  transform an input genomic region into one suitable for the profile.

  :param region:         input region to transform.
  :param window_center:  which part of the input region to center on.
  :param window_size:    how large the resultant region should be.
  :return: a new genomic interval on the same chromosome, centered on the
           <window_center> (e.g. 3' end) of the input region and resized to
           be window_size long.
  """
  if window_center == START:
    return center_start(region)
  elif window_center == END:
    return center_end(region)
  elif window_center == CENTRE:
    return center_middle(region, window_size)
  elif window_center == THREE_PRIME:
    if region.isPositiveStrand():
      return center_end(region, window_size)
    elif region.isNegativeStrand():
      return center_start(region, window_size)
    else:
      raise ValueError("Can't center on three-prime end, no strand.")
  elif window_center == FIVE_PRIME:
    if region.isPositiveStrand():
      return center_start(region, window_size)
    elif region.isNegativeStrand():
      return center_end(region, window_size)
    else:
      raise ValueError("Can't center on five-prime end, no strand.")
  raise ValueError("invalid window center")


###############################################################################
#             HELPER FUNCTIONS FOR BUILDING CONSERVATION PROFILES             #
###############################################################################

def pid(col, ignore_gaps=False):
  """Compute the percent identity of a an alignment column"""
  pass


def conservtion_profile_pid(region, genome_alignment):
  """
  build a conservation profile for the given region using the genome alignment.

  The scores in the profile will be the percent of bases identical to the
  reference sequence.

  :return: a list of the same length as the region where each entry is the
           score of conservation at the corresponding
  """
  res = []
  while len(res) < len(region):
    res.append(0)

  for i in range(region.start, region.end):
    try:
      col = genome_alignment.get_column(region.chrom, i)
      res[i] = pid(col)
    except NoSuchAlignmentColumnError:
      res[i] = 0  # TODO: NaN
    except NoUniqueColumnError:
      res[i] = 0  # TODO: NaN

  return res


def merge_profile(mean_profile, new_profile):
  """
  """
  pass


###############################################################################
#                             MAIN SCRIPT LOGIC                               #
###############################################################################

def processBED(fh, genome_alig, window_size, window_centre, verbose=False):
  """
  Process BED file, produce profile of conservation using whole genome alig.

  :param fh:
  :param genome_alig:   the whole-genome alignment to use to compute
                        conservation scores
  :param window_size:   length of the profile.
  :param window_center: which part of each interval to place at the center
                        of the profile. Acceptable values are in the module
                        constant WINDOW_CENTRE_OPTIONS.
  :param verbose:       if True, output progress messages to stderr.

  :return:
  """
  mean_profile = []
  for e in BEDIterator(fh, verbose=verbose, scoreType=float):
    # figure out which interval to look at...
    region = transform_locus(e)
    new_profile = conservtion_profile_pid(region, genome_alig)
    merge_profile(mean_profile, new_profile)
  return mean_profile


###############################################################################
#                               USER INTERFACE                                #
###############################################################################

def getUI(args):
  programName = os.path.basename(sys.argv[0])
  longDescription = "Given a set of BED intervals, compute a profile of " +\
                    "conservation by averaging over all intervals using a " +\
                    "whole genome alignment to a set of relevent species"
  shortDescription = longDescription

  ui = CLI(programName, shortDescription, longDescription)
  ui.minArgs = 0
  ui.maxArgs = 1
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output to given file, else stdout",
                      required=False, type=str))
  ui.addOption(Option(short="m", long="mafdir", argName="filename",
                      description="directory which contains MAF files",
                      required=True, type=str))
  ui.addOption(Option(short="s", long="species", argName="filename",
                      description="file containing list of assemblies to " +
                                  "use, one per line. All used if omitted",
                      required=False, type=str))
  ui.addOption(Option(short="w", long="window", argName="size",
                      description="size of window to compute around each " +
                                  "interval; " +
                                  str(DEFAULT_WINDOW_SIZE) +
                                  " to use whole interval. " +
                                  "Default " + str(DEFAULT_WINDOW_SIZE),
                      required=False, type=int))
  ui.addOption(Option(short="c", long="centre", argName="location",
                      description="centre window at " + FIVE_PRIME + ", " +
                                  THREE_PRIME + " or " + CENTRE + " " +
                                  "of interval. Ignored if window size " +
                                  "uses the full interval. Default " +
                                  DEFAULT_WINDOW_CENTRE,
                      required=False, type=str))
  ui.addOption(Option(short="t", long="type", argName="str",
                      description="type of scores to compute. Options are:",
                      required=False, type=str))
  ui.addOption(Option(short="v", long="verbose",
                      description="output additional messages to stderr " +
                                  "about run (default: " +
                                  str(DEFAULT_VERBOSITY) + ")",
                      default=DEFAULT_VERBOSITY, required=False))
  ui.addOption(Option(short="h", long="help",
                      description="show this help message ", special=True))
  ui.addOption(Option(short="u", long="test",
                      description="run unit tests ", special=True))

  ui.parseCommandLine(sys.argv[1:])
  return ui


###############################################################################
#                     COMMAND LINE PROCESSING AND DISPATCH                    #
###############################################################################

def main(args):
  # get options and arguments
  ui = getUI()

  # just run unit tests
  if ui.optionIsSet("test"):
    unittest.main(argv=[sys.argv[0]])
    sys.exit()

  # just show help
  if ui.optionIsSet("help"):
    ui.usage()
    sys.exit()
  verbose = (ui.optionIsSet("verbose") is True) or DEFAULT_VERBOSITY

  # make MafDir object
  # mafdir = MafDir(ui.getValue("mafdir"), "maf", "idx")

  # get output handle
  out_fh = sys.stdout
  if ui.optionIsSet("output"):
    out_fh = open(ui.getValue("output"), "w")

  # get input handle
  infh = sys.stdin
  if ui.hasArgument(0):
    infh = open(ui.getArgument(0))

  # get window size...
  windowSize = DEFAULT_WINDOW_SIZE
  if ui.optionIsSet("window"):
    windowSize = ui.getValue("window")

  # get the window anchoring location
  windowCentre = DEFAULT_WINDOW_CENTRE
  if ui.optionIsSet("centre"):
    windowCentre = ui.getValue("centre")
    if windowCentre not in WINDOW_CENTRE_OPTIONS:
      sys.stderr.write("un-recognised window anchor position: " +
                       str(windowCentre) + "\n")
      sys.exit(1)

  # get set of species to use
  species = None
  if ui.optionIsSet("species"):
    species = [l.strip() for l in open(ui.getValue("species"))
               if l.strip() != ""]

  # processBED(infh, out_fh, mafdir, windowSize, windowCentre, species, ref,
  #           verbose)


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    main()
