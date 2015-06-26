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
import StringIO

# pyokit imports
from pyokit.interface.cli import CLI
from pyokit.interface.cli import Option
from pyokit.io.bedIterators import BEDIterator
from pyokit.datastruct.genomeAlignment import NoSuchAlignmentColumnError
from pyokit.datastruct.genomeAlignment import NoUniqueColumnError
from pyokit.datastruct import sequence
from pyokit.statistics.online import RollingMean
from pyokit.io.genomeAlignment import genome_alignment_iterator
from pyokit.datastruct.genomeAlignment import GenomeAlignment


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
  Center a region on its end and expand it to window_size bases.

  :return: the new region.
  """
  res = copy.copy(r)
  res.start = res.end - window_size / 2
  res.end = res.start + window_size
  return res


def center_middle(r, window_size):
  """
  Center a region on its middle and expand it to window_size bases.

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
  if window_center == CENTRE:
    region.transform_center(window_size)
  else:
    raise ValueError("Don't know how to do this transformation: " +
                     window_center)


###############################################################################
#             HELPER FUNCTIONS FOR BUILDING CONSERVATION PROFILES             #
###############################################################################

def pid(col, ignore_gaps=False):
  """
    Compute the percent identity of a an alignment column.

    Define PID as the frequency of the most frequent nucleotide in the column.

    :param col:         an alignment column; a dictionary where keys are seq.
                        names and values are the nucleotide in the column for
                        that sequence.
    :param ignore_gaps: if True, do not count gaps towards the total number of
                        sequences in the column (i.e. the denominator of the
                        fraction).
    :raise ValueError: if the column contains only gaps.
  """
  hist = {}
  total = 0
  found_non_gap = False
  for v in col.values():
    if v == sequence.GAP_CHAR:
      if ignore_gaps:
        continue
      else:
        total += 1
    else:
      found_non_gap = True
      if v not in hist:
        hist[v] = 0
      hist[v] += 1
      total += 1
  if not found_non_gap:
    raise ValueError("Cannot determine PID of column with only gaps")
  return max(hist.values()) / float(total)


def conservtion_profile_pid(region, genome_alignment):
  """
  build a conservation profile for the given region using the genome alignment.

  The scores in the profile will be the percent of bases identical to the
  reference sequence.

  :return: a list of the same length as the region where each entry is the
           PID at the corresponding locus.
  """
  res = []
  s = region.start if region.isPositiveStrand() else region.end - 1
  e = region.end if region.isPositiveStrand() else region.start - 1
  step = 1 if region.isPositiveStrand() else -1
  for i in range(s, e, step):
    try:
      col = genome_alignment.get_column(region.chrom, i)
      res.append(pid(col))
    except NoSuchAlignmentColumnError:
      res.append(None)
    except NoUniqueColumnError:
      res.append(None)

  return res


def merge_profile(mean_profile, new_profile):
  """Add a new list of values to a list of rolling means."""
  for i in range(0, len(mean_profile)):
    if new_profile[i] is None:
      continue
    mean_profile[i].add(new_profile[i])


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
  while len(mean_profile) < window_size:
    mean_profile.append(RollingMean())

  for e in BEDIterator(fh, verbose=verbose, scoreType=float):
    # figure out which interval to look at...
    transform_locus(e, window_centre, window_size)
    new_profile = conservtion_profile_pid(e, genome_alig)
    merge_profile(mean_profile, new_profile)
  return [m.mean for m in mean_profile]


###############################################################################
#                               USER INTERFACE                                #
###############################################################################

def getUI(args):
  """Build and return user interface object for this script."""
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
  """Process the command line arguments of this script and dispatch."""
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

  # processBED(infh, out_fh, mafdir, windowSize, windowCentre, species, ref,
  #           verbose)


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestConservationProfile(unittest.TestCase):

  """Unit tests for this script."""

  def setUp(self):
    """Set up a genome alignment and set of query regions for testing."""
    b1_A_seq = "CTGATGCAGTC-"
    b1_B_seq = "TA-ATGCA-ATG"
    b1_B_qul = "gg-ggggg-ggg"
    b1_C_seq = "TGGT-GCAGTAA"
    b1_C_qul = "gg-ggggg-ggg"
    b1 = "a score=28452.00\n" +\
         "s A.chr1 10 11 + 51323347 " + b1_A_seq + "\n" +\
         "s B.chr1 50 10 + 51313456 " + b1_B_seq + "\n" +\
         "q B.chr1                  " + b1_B_qul + "\n" +\
         "s C.chr2 31 11 + 50314486 " + b1_C_seq + "\n" +\
         "q C.chr2                  " + b1_C_qul
    b2_A_seq = "GTAGTAGC-"
    b2_B_seq = "-CAAT-GCA"
    b2_B_qul = "-gggg-ggg"
    b2_C_seq = "-A-TAAGCC"
    b2_C_qul = "-g-gggggg"
    b2 = "a score=45845.2456\n" +\
         "s A.chr1 88 8 + 51323347 " + b2_A_seq + "\n" +\
         "s B.chr1 7  7 + 51313456 " + b2_B_seq + "\n" +\
         "q B.chr1                 " + b2_B_qul + "\n" +\
         "s C.chr1 62 7 + 50314486 " + b2_C_seq + "\n" +\
         "q C.chr1                 " + b2_C_qul
    self.maf1 = b1 + "\n\n" + b2

    b3_A_seq = "G-CCGATGC"
    b3_B_seq = "ACCC-CTGA"
    b3_B_qul = "gggg-gggg"
    b3_C_seq = "ACCC-GGGA"
    b3_C_qul = "gggg-gggg"
    b3 = "a score=15489.458\n" +\
         "s A.chrX 10 8 + 53347 " + b3_A_seq + "\n" +\
         "s B.chrX 13 8 + 53456 " + b3_B_seq + "\n" +\
         "q B.chrX              " + b3_B_qul + "\n" +\
         "s C.chrX 9  8 - 50486 " + b3_C_seq + "\n" +\
         "q C.chrX              " + b3_C_qul
    b4_A_seq = "GGAGTTA"
    b4_B_seq = "-GA-TTA"
    b4_B_qul = "-gg-ggg"
    b4 = "a score=14489.458\n" +\
         "s A.chr7 10 7 + 53347 " + b4_A_seq + "\n" +\
         "s B.chr7 33 5 + 53456 " + b4_B_seq + "\n" +\
         "q B.chr7              " + b4_B_qul + "\n" +\
         "e C.chr7 54  7 + 50486 I"
    self.maf2 = b3 + "\n\n" + b4

    self.roi = "chr1" + "\t" + "14" + "\t" + "18" + "\t" + "X" + "\t" + "0" +\
               "\t" + "+\n" +\
               "chr1" + "\t" + "19" + "\t" + "22" + "\t" + "X" + "\t" + "0" +\
               "\t" + "+\n" +\
               "chr1" + "\t" + "88" + "\t" + "92" + "\t" + "X" + "\t" + "0" +\
               "\t" + "+\n" +\
               "chr1" + "\t" + "94" + "\t" + "98" + "\t" + "X" + "\t" + "0" +\
               "\t" + "-\n" +\
               "chrX" + "\t" + "11" + "\t" + "14" + "\t" + "X" + "\t" + "0" +\
               "\t" + "-\n" +\
               "chrX" + "\t" + "15" + "\t" + "18" + "\t" + "X" + "\t" + "0" +\
               "\t" + "+\n" +\
               "chr7" + "\t" + "10" + "\t" + "12" + "\t" + "X" + "\t" + "0" +\
               "\t" + "+\n"

  def test_conservation_profile_pid(self):
    """Test getting conservation profiles (PID) from genome alignments"""
    m = StringIO.StringIO(self.maf1 + "\n\n" + self.maf2)
    ga = GenomeAlignment([x for x in genome_alignment_iterator(m, "A")])

    # here we're testing directly; the actual script will adjust the regions
    # before it does this step.
    expect_raw = [[0.66666666, 1.00000000, 1.00000000, 1.0000000],
                  [0.66666666, 0.33333333, None],
                  [0.33333333, 0.33333333, 0.66666666, 0.33333333],
                  [None, None, 1.00000000, 1.00000000],
                  [0.33333333, 1.00000000, 1.00000000],
                  [0.66666666, 1.00000000, 0.66666666],
                  [0.33333333, 0.66666666]]
    in_regions = [r for r in BEDIterator(StringIO.StringIO(self.roi))]
    res = [conservtion_profile_pid(r, ga) for r in in_regions]
    self.assertEqual(len(expect_raw), len(res))
    for i in range(0, len(expect_raw)):
      self.assertEqual(len(expect_raw[i]), len(res[i]))
      for j in range(0, len(expect_raw[i])):
        self.assertAlmostEqual(expect_raw[i][j], res[i][j])

    # now we test with the adjusted regions
    for l in in_regions:
      transform_locus(l, CENTRE, 4)
    res_adjusted = [conservtion_profile_pid(r, ga) for r in in_regions]
    expect_adjusted = [[0.66666666, 1.00000000, 1.00000000, 1.0000000],
                       [0.66666666, 0.33333333, None, None],
                       [0.33333333, 0.33333333, 0.66666666, 0.33333333],
                       [None, None, 1.00000000, 1.00000000],
                       [0.33333333, 1.00000000, 1.00000000, 0.66666666],
                       [0.66666666, 1.00000000, 0.66666666, None],
                       [None, 0.33333333, 0.66666666, 0.66666666]]
    self.assertEqual(len(expect_adjusted), len(res_adjusted))
    for i in range(0, len(expect_adjusted)):
      self.assertEqual(len(expect_adjusted[i]), len(res_adjusted[i]))
      for j in range(0, len(expect_adjusted[i])):
        self.assertAlmostEqual(expect_adjusted[i][j], res_adjusted[i][j])

  def test_process(self):
    """Bypass the UI and test the main logic function directly."""
    m = StringIO.StringIO(self.maf1 + "\n\n" + self.maf2)
    ga = GenomeAlignment([x for x in genome_alignment_iterator(m, "A")])
    prof = processBED(StringIO.StringIO(self.roi), ga, 4, CENTRE)
    expect = [0.53333333, 0.66666666, 0.83333333, 0.73333333]
    self.assertEqual(len(prof), len(expect))
    for i in range(0, len(expect)):
      self.assertAlmostEqual(prof[i], expect[i])


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    main()
