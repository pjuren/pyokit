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
import functools

# for testing
import mock

# pyokit imports -- IO
from pyokit.io.genomeAlignment import genome_alignment_iterator
from pyokit.io.genomeAlignment import build_genome_alignment_from_file
from pyokit.io.genomeAlignment import load_just_in_time_genome_alignment
from pyokit.io.indexedFile import IndexedFile
from pyokit.io.bedIterators import BEDIterator
from pyokit.io.bedIterators import ITERATOR_SORTED_START
# pyokit improts - UI
from pyokit.interface.cli import CLI
from pyokit.interface.cli import Option
# pyokit imports data structures
from pyokit.datastruct.genomeAlignment import NoSuchAlignmentColumnError
from pyokit.datastruct.genomeAlignment import NoUniqueColumnError
from pyokit.datastruct.genomeAlignment import JustInTimeGenomeAlignmentBlock
from pyokit.datastruct import sequence
from pyokit.datastruct.genomeAlignment import GenomeAlignment
from pyokit.datastruct.multipleAlignment import MissingSequenceHandler
# pyokit imports -- statistics
from pyokit.statistics.online import RollingMean


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


def conservtion_profile_pid(region, genome_alignment,
                            mi_seqs=MissingSequenceHandler.TREAT_AS_ALL_GAPS,
                            species=None):
  """
  build a conservation profile for the given region using the genome alignment.

  The scores in the profile will be the percent of bases identical to the
  reference sequence.

  :param miss_seqs: how to treat sequence with no actual sequence data for
                    the column.
  :return: a list of the same length as the region where each entry is the
           PID at the corresponding locus.
  """
  res = []
  s = region.start if region.isPositiveStrand() else region.end - 1
  e = region.end if region.isPositiveStrand() else region.start - 1
  step = 1 if region.isPositiveStrand() else -1
  for i in range(s, e, step):
    try:
      col = genome_alignment.get_column(region.chrom, i, mi_seqs, species)
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

def processBED(fh, genome_alig, window_size, window_centre,
               mi_seqs=MissingSequenceHandler.TREAT_AS_ALL_GAPS, species=None,
               verbose=False):
  """
  Process BED file, produce profile of conservation using whole genome alig.

  :param fh:
  :param genome_alig:   the whole-genome alignment to use to compute
                        conservation scores
  :param window_size:   length of the profile.
  :param window_center: which part of each interval to place at the center
                        of the profile. Acceptable values are in the module
                        constant WINDOW_CENTRE_OPTIONS.
  :param miss_seqs:     how to treat sequence with no actual sequence data for
                        the column.
  :param verbose:       if True, output progress messages to stderr.

  :return:
  """
  mean_profile = []
  while len(mean_profile) < window_size:
    mean_profile.append(RollingMean())

  for e in BEDIterator(fh, verbose=verbose, scoreType=float,
                       sortedby=ITERATOR_SORTED_START):
    # figure out which interval to look at...
    transform_locus(e, window_centre, window_size)
    new_profile = conservtion_profile_pid(e, genome_alig, mi_seqs, species)
    merge_profile(mean_profile, new_profile)
  return [m.mean for m in mean_profile]


###############################################################################
#                               USER INTERFACE                                #
###############################################################################

def getUI(prog_name, args):
  """Build and return user interface object for this script."""
  longDescription = "Given a set of BED intervals, compute a profile of " +\
                    "conservation by averaging over all intervals using a " +\
                    "whole genome alignment to a set of relevent species." +\
                    "\n\n" +\
                    "Usage: " + prog_name + " [options] regions.bed " +\
                    "genome-alig species" +\
                    "\n\n" +\
                    "genome-alig can be either a single MAF file, or a " +\
                    "directory of MAF files. In the latter case, the " +\
                    "directory may also optionally contain index files for " +\
                    "the alignment files."
  shortDescription = longDescription

  ui = CLI(prog_name, shortDescription, longDescription)
  # gotta have two args -- MAF dir/file and BED regions.
  # Input by stdin not allowed
  ui.minArgs = 3
  ui.maxArgs = 4
  ui.addOption(Option(short="o", long="output", argName="filename",
                      description="output to given file, else stdout",
                      required=False, type=str))
  ui.addOption(Option(short="w", long="window", argName="size",
                      description="size of window to compute around each " +
                                  "interval; " +
                                  str(DEFAULT_WINDOW_SIZE) +
                                  " to use whole interval. " +
                                  "Default " + str(DEFAULT_WINDOW_SIZE),
                      required=False, type=int))
  ui.addOption(Option(short="e", long="extensions", argName="extension",
                      description="if genome-alig specifies a directory, " +
                                  "treat files with this extension as " +
                                  "alignment files.", required=False,
                                  type=str))
  ui.addOption(Option(short="i", long="index-extensions", argName="extension",
                      description="if genome-alig specifies a directory, " +
                                  "treat files with this extension as " +
                                  "index files for alignments.",
                                  required=False, type=str))
  ui.addOption(Option(short="f", long="fail-no-index",
                      description="fail if an alignment file without an " +
                                  "index is found; otherwise index-less " +
                                  "alignment files are loaded whole (which " +
                                  "might be slow if they're large, and " +
                                  "might require a lot of memory)",
                      default=False, required=False))
  ui.addOption(Option(short="m", long="missing", argName="strategy",
                      description="how to treat missing sequences in " +
                                  "blocks. Options are " +
                                  ", ".join([str(x.name) for x in
                                             MissingSequenceHandler]),
                                  required=False, type=str))
  ui.addOption(Option(short="s", long="species", argName="species",
                      description="consider only these species. Default is " +
                                  "all.", required=False, type=str))
  ui.addOption(Option(short="v", long="verbose",
                      description="output additional messages to stderr " +
                                  "about run (default: " +
                                  str(DEFAULT_VERBOSITY) + ")",
                      default=DEFAULT_VERBOSITY, required=False))
  ui.addOption(Option(short="h", long="help",
                      description="show this help message ", special=True))
  ui.addOption(Option(short="u", long="test",
                      description="run unit tests ", special=True))

  ui.parseCommandLine(args)
  return ui


###############################################################################
#                     COMMAND LINE PROCESSING AND DISPATCH                    #
###############################################################################

def main(args, prog_name):
  """Process the command line arguments of this script and dispatch."""
  # get options and arguments
  ui = getUI(prog_name, args)

  # just run unit tests
  if ui.optionIsSet("test"):
    unittest.main(argv=[sys.argv[0]])
    sys.exit()

  # just show help
  if ui.optionIsSet("help"):
    ui.usage()
    sys.exit()
  verbose = (ui.optionIsSet("verbose") is True) or DEFAULT_VERBOSITY

  # get output handle
  out_fh = sys.stdout
  if ui.optionIsSet("output"):
    out_fh = open(ui.getValue("output"), "w")

  # get window size...
  window_size = DEFAULT_WINDOW_SIZE
  if ui.optionIsSet("window"):
    window_size = ui.getValue("window")

  # get the window anchoring location
  windowCentre = DEFAULT_WINDOW_CENTRE
  if ui.optionIsSet("centre"):
    windowCentre = ui.getValue("centre")
    if windowCentre not in WINDOW_CENTRE_OPTIONS:
      sys.stderr.write("un-recognised window anchor position: " +
                       str(windowCentre) + "\n")
      sys.exit(1)

  args = ui.getAllArguments()
  assert(len(args) == 3 or len(args) == 4)
  region_fn = ui.getArgument(0)
  ga_path = ui.getArgument(1)
  index_fn = None
  if len(args) == 3:
    spec = ui.getArgument(2)
  else:
    index_fn = ui.getArgument(2)
    spec = ui.getArgument(3)

  extensions = (ui.getValue("extensions").strip().split(",") if
                ui.optionIsSet("extensions") else None)
  index_extensions = (ui.getValue("index-extensions").strip().split(",") if
                      ui.optionIsSet("index-extensions") else None)
  fail_no_index = ui.optionIsSet("fail-no-index")
  mi_seqs = (MissingSequenceHandler[ui.getValue("missing")]
             if ui.optionIsSet("missing")
             else MissingSequenceHandler.TREAT_AS_ALL_GAPS)
  species = ([x.strip() for x in ui.getValue("species").split(",")] if
             ui.optionIsSet("species") else None)

  # build the genome alignment
  alig = (load_just_in_time_genome_alignment(ga_path, spec, extensions,
                                             index_extensions, fail_no_index,
                                             verbose)
          if os.path.isdir(ga_path)
          else build_genome_alignment_from_file(ga_path, spec, index_fn,
                                                verbose))

  # get the profile and write it to the output stream
  profile = processBED(open(region_fn), alig, window_size, CENTRE,
                       mi_seqs, species, verbose)
  out_fh.write("\n\n" + ", ".join(str(x) for x in profile) + "\n")


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

def _build_index(in_strng, ref_spec):
  idx_strm = StringIO.StringIO()
  bound_iter = functools.partial(genome_alignment_iterator,
                                 reference_species=ref_spec)
  hash_func = JustInTimeGenomeAlignmentBlock.build_hash
  idx = IndexedFile(StringIO.StringIO(in_strng), bound_iter, hash_func)
  idx.write_index(idx_strm)
  idx_strm.seek(0)  # seek to the start
  return idx_strm


def _build_open_side_effect(lookup):
  def open_side_effect(*args, **kwargs):
    if not isinstance(args[0], basestring):
      raise TypeError()
    fn = args[0].strip()
    if fn in lookup:
      return lookup[fn]
    raise IOError("No such file: " + args[0])
  return open_side_effect


def _build_isfile_side_effect(files_lookup, dirs_lookup):
  def isfile_side_effect(*args, **kwargs):
    fn = args[0].strip()
    if fn in files_lookup:
      return True
    if fn in dirs_lookup:
      return False
    raise IOError("No such file or directory: " + fn)
  return isfile_side_effect


def _build_isdir_side_effect(files_lookup, dirs_lookup):
  def isdir_side_effect(*args, **kwargs):
    fn = args[0].strip()
    if fn in files_lookup:
      return False
    if fn in dirs_lookup:
      return True
    raise IOError("No such file or directory: " + fn)
  return isdir_side_effect


class ConsProfileTestHelper(object):

  """A static class containing a set of genome alignment blocks for testing."""

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
  maf1 = b1 + "\n\n" + b2

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
  maf2 = b3 + "\n\n" + b4

  roi = "\t".join(["chr1", "14", "18", "X", "0", "+\n"]) +\
        "\t".join(["chr1", "19", "22", "X", "0", "+\n"]) +\
        "\t".join(["chr1", "88", "92", "X", "0", "+\n"]) +\
        "\t".join(["chr1", "94", "98", "X", "0", "-\n"]) +\
        "\t".join(["chr7", "10", "12", "X", "0", "+\n"]) +\
        "\t".join(["chrX", "11", "14", "X", "0", "-\n"]) +\
        "\t".join(["chrX", "15", "18", "X", "0", "+\n"])


class TestConservationProfileIndvFiles(unittest.TestCase):

  """Unit tests for this script."""

  def setUp(self):
    """Set up a genome alignment and set of query regions for testing."""
    self.b1 = ConsProfileTestHelper.b1
    self.b2 = ConsProfileTestHelper.b2
    self.maf1 = ConsProfileTestHelper.maf1
    self.b3 = ConsProfileTestHelper.b3
    self.b4 = ConsProfileTestHelper.b4
    self.maf2 = ConsProfileTestHelper.maf2
    self.roi = ConsProfileTestHelper.roi

  def test_conservation_profile_pid(self):
    """Test getting conservation profiles (PID) from genome alignments."""
    m = StringIO.StringIO(self.maf1 + "\n\n" + self.maf2)
    ga = GenomeAlignment([x for x in genome_alignment_iterator(m, "A")])

    # here we're testing directly; the actual script will adjust the regions
    # before it does this step.
    expect_raw = [[0.66666666, 1.00000000, 1.00000000, 1.0000000],
                  [0.66666666, 0.33333333, None],
                  [0.33333333, 0.33333333, 0.66666666, 0.33333333],
                  [None, None, 1.00000000, 1.00000000],
                  [0.33333333, 0.66666666],
                  [0.33333333, 1.00000000, 1.00000000],
                  [0.66666666, 1.00000000, 0.66666666]]
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
                       [None, 0.33333333, 0.66666666, 0.66666666],
                       [0.33333333, 1.00000000, 1.00000000, 0.66666666],
                       [0.66666666, 1.00000000, 0.66666666, None]]
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

  @mock.patch('__builtin__.open')
  def test_full_with_UI(self, mock_open):
    """Test the full script, including the UI."""
    output_stream = StringIO.StringIO()

    def open_side_effect(*args, **kwargs):
      if args[0] == "in.maf":
        return StringIO.StringIO(self.maf1 + "\n\n" + self.maf2)
      if args[0] == "in.bed":
        return StringIO.StringIO(self.roi)
      if args[0] == "out.txt":
        return output_stream
      raise IOError("No such file")

    mock_open.side_effect = open_side_effect

    main(["-w", "4", "-o", "out.txt", "in.bed", "in.maf", "A"], "cons_profile")
    vals = [float(x) for x in output_stream.getvalue().split(",")]
    expect = [0.53333333, 0.66666666, 0.83333333, 0.73333333]
    self.assertEqual(len(expect), len(vals))
    for i in range(0, len(vals)):
      self.assertAlmostEqual(vals[i], expect[i])


class TestConservationProfileDirectory(unittest.TestCase):

  """Test using a directory of files for alignments."""

  def setUp(self):
    """Prepare a directory structure and files for tests."""
    self.b1 = ConsProfileTestHelper.b1
    self.b2 = ConsProfileTestHelper.b2
    self.b3 = ConsProfileTestHelper.b3
    self.b4 = ConsProfileTestHelper.b4
    self.roi = ConsProfileTestHelper.roi

    # set up a directory structure
    jn = os.path.join
    self.lk_up = {jn("the_dir", "chr1:10-21.maf"): StringIO.StringIO(self.b1),
                  jn("the_dir", "chr1:88-96.maf"): StringIO.StringIO(self.b2),
                  jn("the_dir", "chrX:10-17.maf"): StringIO.StringIO(self.b3),
                  jn("the_dir", "chr7:10-17.maf"): StringIO.StringIO(self.b4),
                  jn("the_dir", "chr1:10-21.idx"): _build_index(self.b1, "A"),
                  jn("the_dir", "chr1:88-96.idx"): _build_index(self.b2, "A"),
                  jn("the_dir", "chrX:10-17.idx"): _build_index(self.b3, "A"),
                  jn("the_dir", "chr7:10-17.idx"): _build_index(self.b4, "A"),
                  "in.bed": StringIO.StringIO(self.roi)}

  @mock.patch('os.path.isdir')
  @mock.patch('os.listdir')
  @mock.patch('os.path.isfile')
  @mock.patch('__builtin__.open')
  def test_full_UI_just_in_time_genome_alig(self, mock_open, mock_isfile,
                                            mock_listdir, mock_isdir):
    """Test using a JIT genome alignment built from a directory."""
    # b1 --> chr1 10 21  ;  b2 --> chr1 88 96
    # b3 --> chrX 10 17  ;  b4 --> chr7 10 17
    mock_listdir.return_value = ["chr1:10-21.maf", "chr1:88-96.maf",
                                 "chrX:10-17.maf", "chr7:10-17.maf",
                                 "chr1:10-21.idx", "chr1:88-96.idx",
                                 "chrX:10-17.idx", "chr7:10-17.idx",
                                 "some_sub_dir"]
    output_stream = StringIO.StringIO()
    dir_lkup = set(["the_dir", os.path.join("the_dir", "some_sub_dir")])
    file_lkup = self.lk_up
    output_stream = StringIO.StringIO()
    file_lkup["out.txt"] = output_stream

    mock_open.side_effect = _build_open_side_effect(file_lkup)
    mock_isfile.side_effect = _build_isfile_side_effect(file_lkup, dir_lkup)
    mock_isdir.side_effect = _build_isdir_side_effect(file_lkup, dir_lkup)

    main(["-w", "4", "-o", "out.txt", "-e", ".maf", "-i", ".idx", "in.bed",
          "the_dir", "A"], "cons_profile")
    vals = [float(x) for x in output_stream.getvalue().split(",")]
    expect = [0.53333333, 0.66666666, 0.83333333, 0.73333333]
    self.assertEqual(len(expect), len(vals))
    for i in range(0, len(vals)):
      self.assertAlmostEqual(vals[i], expect[i])


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    main(sys.argv[1:], sys.argv[0])
