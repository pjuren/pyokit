"""
Date of Creation: 11th Dec 2014.

Description:   Code for reading/writing genome alignments.

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

# standard python imports
import StringIO
import os
import os.path
import functools

# import for unit testing
import mock
import unittest

# pyokit imports
from pyokit.datastruct.genomeAlignment import GenomeAlignmentBlock
from pyokit.datastruct.genomeAlignment import JustInTimeGenomeAlignmentBlock
from pyokit.datastruct.genomeAlignment import JustInTimeGenomeAlignment
from pyokit.datastruct.genomeAlignment import GenomeAlignment
from pyokit.datastruct.sequence import Sequence
from pyokit.datastruct.sequence import UnknownSequence
from pyokit.io import maf
from pyokit.io.indexedFile import IndexedFile
from pyokit.util.progressIndicator import ProgressIndicator
from pyokit.io.ioError import PyokitIOError


# def genome_alignment_block_hash(b):
#  """Hash a genome alignment block to unique ID: its location in the genome."""
#  return b.chrom + "\t" + str(b.start) + "\t" + str(b.end)

###############################################################################
#                              HELPER FUNCTIONS                               #
###############################################################################

def __trim_extension_dot(ext):
  """trim leading dots from extension."""
  if ext is None:
    return None
  if ext == "":
    return ""
  while ext[0] == ".":
    ext = ext[1:]
  return ext


def __trim_extensions_dot(exts):
  """trim leading dots from extensions and drop any empty strings."""
  if exts is None:
    return None
  res = []
  for i in range(0, len(exts)):
    if exts[i] == "":
      continue
    res.append(__trim_extension_dot(exts[i]))
  return res


def __split_genomic_interval_filename(fn):
  """
  Split a filename of the format chrom:start-end.ext or chrom.ext (full chrom).

  :return: tuple of (chrom, start, end) -- 'start' and 'end' are None if not
           present in the filename.
  """
  if fn is None or fn == "":
    raise ValueError("invalid filename: " + str(fn))
  fn = ".".join(fn.split(".")[:-1])
  parts = fn.split(":")
  if len(parts) == 1:
    return (parts[0].strip(), None, None)
  else:
    r_parts = parts[1].split("-")
    if len(r_parts) != 2:
      raise ValueError("Invalid filename: " + str(fn))
    return (parts[0].strip(), int(r_parts[0]), int(r_parts[1]))


###############################################################################
#                MANAGING AN ON-DISK MULTI-FILE GENOME ALIGNMENT              #
###############################################################################

def load_just_in_time_genome_alignment(path, ref_spec, extensions=None,
                                       index_exts=None, fail_no_index=True,
                                       verbose=False):
    """Load a just-in-time genome alignment from a directory."""
    extensions = __trim_extensions_dot(extensions)
    index_exts = __trim_extensions_dot(index_exts)

    partial_chrom_files = {}
    whole_chrom_files = {}
    for fn in os.listdir(path):
      pth = os.path.join(path, fn)
      if os.path.isfile(pth):
        base, ext = os.path.splitext(pth)
        ext = __trim_extension_dot(ext)
        if extensions is None or ext in extensions:
          idx_path = __find_index(pth, index_exts)
          if idx_path is None and fail_no_index:
            raise PyokitIOError("No index file for " + fn)
          chrom, start, end = __split_genomic_interval_filename(fn)
          assert((start is None and end is None) or
                 (start is not None and end is not None))
          if start is None:
            if chrom in whole_chrom_files:
              raise PyokitIOError("multiple files for chrom " + chrom)
            whole_chrom_files[chrom] = (pth, idx_path)
          else:
            k = (chrom, start, end)
            if k in partial_chrom_files:
              raise PyokitIOError("multiple files for " + str(k))
            partial_chrom_files[k] = (pth, idx_path)

    def factory(k):
      pth, idx = k
      return build_genome_alignment_from_file(pth, ref_spec, idx, verbose)

    return JustInTimeGenomeAlignment(whole_chrom_files, partial_chrom_files,
                                     factory)


###############################################################################
#                   BULK LOADING GENOME ALIGNMENTS FROM DISK                  #
###############################################################################

def __find_index(alig_file_pth, idx_extensions):
  """
  Find an index file for a genome alignment file in the same directory.

  :param alig_file_path: path to the alignment file.
  :param idx_extensions: check for index files with these extensions
  :return: path to first index file that matches the name of the alignment file
           and has one of the specified extensions.
  """
  if idx_extensions is None:
    return None
  base, ext = os.path.splitext(alig_file_pth)
  for idx_ext in idx_extensions:
    candidate = base + os.extsep + idx_ext
    if os.path.isfile(candidate):
      return candidate
  return None


def build_genome_alignment_from_directory(d_name, ref_spec, extensions=None,
                                          index_exts=None, fail_no_index=False,
                                          verbose=False):
  """
  build a genome aligment by loading all files in a directory.

  Fiel without indexes are loaded immediately; those with indexes are
  loaded on-demand. Not recursive (i.e. subdirectories are not parsed).

  :param d_name:        directory to load from.
  :param ref_spec:      which species in the alignemnt files is the reference?
  :param extensions:    list or set of acceptable extensions; treat any files
                        with these extensions as part of the alignment. If None,
                        treat any file which has an extension that is NOT in
                        index_extensions as part of the alignment.
  :param index_exts:    treat any files with these extensions as index files.
  :param fail_no_index: fail if index extensions are provided and an alignment
                        file has no index file.
  """
  if index_exts is None and fail_no_index:
    raise ValueError("Failure on no index specified for loading genome " +
                     "alignment, but no index extensions specified")

  blocks = []
  for fn in os.listdir(d_name):
    pth = os.path.join(d_name, fn)
    if os.path.isfile(pth):
      base, ext = os.path.splitext(pth)
      if extensions is None or ext in extensions:
        idx_path = __find_index(pth, index_exts)
        if idx_path is None and fail_no_index:
          raise PyokitIOError("No index file for " + fn)
        for b in genome_alignment_iterator(pth, ref_spec, idx_path):
          blocks.append(b)
  return GenomeAlignment(blocks)


def build_genome_alignment_from_file(ga_path, ref_spec, idx_path=None,
                                     verbose=False):
  """
  build a genome alignment by loading from a single MAF file.

  :param ga_path:  the path to the file to load.
  :param ref_spec: which species in the MAF file is the reference?
  :param idx_path: if provided, use this index to generate a just-in-time
                   genome alignment, instead of loading the file immediately.
  """
  blocks = []
  if (idx_path is not None):
    bound_iter = functools.partial(genome_alignment_iterator,
                                   reference_species=ref_spec)
    hash_func = JustInTimeGenomeAlignmentBlock.build_hash
    factory = IndexedFile(None, bound_iter, hash_func)
    factory.read_index(idx_path, ga_path, verbose=verbose)

    pind = None
    for k in factory:
      if verbose:
        if pind is None:
          total = len(factory)
          pind = ProgressIndicator(totalToDo=total, messagePrefix="completed",
                                   messageSuffix="building alignment blocks ")
        pind.done += 1
        pind.showProgress()
      blocks.append(JustInTimeGenomeAlignmentBlock(factory, k))
  else:
    for b in genome_alignment_iterator(ga_path, ref_spec, verbose=verbose):
      blocks.append(b)
  return GenomeAlignment(blocks)


###############################################################################
#                                ITERATORS                                    #
###############################################################################

def genome_alignment_iterator(fn, reference_species, index_friendly=False,
                              verbose=False):
  """
  build an iterator for an MAF file of genome alignment blocks.

  :param fn:                 filename or stream-like object to iterate over.
  :param reference_species:  which species in the alignment should be treated
                             as the reference?
  :param index_friendly:     if True, buffering is disabled to support using
                             the iterator to build an index.

  :return an iterator that yields GenomeAlignment objects
  """
  kw_args = {"reference_species": reference_species}
  for e in maf.maf_iterator(fn, index_friendly=index_friendly,
                            yield_class=GenomeAlignmentBlock,
                            yield_kw_args=kw_args):
    yield e


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################

class WrappedStringIO(object):

  """A wrapper for a stringIO object that tracks of the number of reads."""

  def __init__(self, in_str):
    """Constructor for the WrappedStringIO class. in_strm = stream to wrap."""
    self.io = StringIO.StringIO(in_str)
    self.reads = 0

  def tell(self):
    """Wrapper for the StringIO.tell function."""
    return self.io.tell()

  def seek(self, x):
    """wrapper for the StringIO.seek function."""
    self.io.seek(x)

  def __iter__(self):
    """Wrapper for the StringIO iteration functionality."""
    for l in self.io:
      self.reads += 1
      yield l
    raise StopIteration()

  def read(self):
    """Wrapper for the StringIO read function."""
    self.reads += 1
    return self.io.read()


def _build_index(maf_strm, ref_spec):
  """Build an index for a MAF genome alig file and return StringIO of it."""
  idx_strm = StringIO.StringIO()
  bound_iter = functools.partial(genome_alignment_iterator,
                                 reference_species=ref_spec)
  hash_func = JustInTimeGenomeAlignmentBlock.build_hash
  idx = IndexedFile(maf_strm, bound_iter, hash_func)
  idx.write_index(idx_strm)
  idx_strm.seek(0)  # seek to the start
  return idx_strm


class GATestHelper(object):

  """Helper for tests involving genome alignments in concrete MAF syntax."""

  # this defines a single genome alignmnet block
  b1_hg19_seq = "atctccaagagggcataaaacac-tgagtaaacagctcttttatatgtgtttcctgga"
  b1_panTro_s = "atctccaagagggcataaaacac-tgagtaaacagctctt--atatgtgtttcctgga"
  b1_panTro_q = "99999999999999999999999-9999999999999999--9999999999999999"
  b1_tarSyr_s = "atctccaagagggctgaaaatgc-caaatga-----------tcacacgtttcctgga"
  b1_tarSyr_q = "79295966999999999999998-9999799-----------9999999999765775"
  b1_tupBel_s = "ttcaggaagggggcccaaaacgcttgagtggtcagctctta-ttttgcgtttactgga"
  b1_tupBel_q = "79648579699867994997775679665662767577569-6998745597677632"
  b1 = "a score=28680.000000\n" +\
       "s hg19.chr22             1711 57 + 51304566 " + b1_hg19_seq + "\n" +\
       "s panTro2.chrUn          1110 59 + 58616431 " + b1_panTro_s + "\n" +\
       "q panTro2.chrUn                             " + b1_panTro_q + "\n" +\
       "i panTro2.chrUn          C 0 C 0                          " + "\n" +\
       "s tarSyr1.scaffold_5923  2859 50 -     8928 " + b1_tarSyr_s + "\n" +\
       "q tarSyr1.scaffold_5923                     " + b1_tarSyr_q + "\n" +\
       "i tarSyr1.scaffold_5923  N 0 C 0                          " + "\n" +\
       "s tupBel1.scaffold_803   33686 61 +   85889 " + b1_tupBel_s + "\n" +\
       "q tupBel1.scaffold_803                      " + b1_tupBel_q + "\n" +\
       "i tupBel1.scaffold_803   I 1 C 0                          " + "\n" +\
       "e mm4.chr6            53310102 58 + 151104725 I"
  b1_hg19 = Sequence("hg19.chr22", b1_hg19_seq, 1711, 1768,
                     "+", 51302798)
  b1_panTro = Sequence("panTro2.chrUn", b1_panTro_s, 1110, 1169, "+",
                       58616431 - 1169,
                       {maf.QUALITY_META_KEY: b1_panTro_q,
                        maf.LEFT_STATUS_KEY: "C",
                        maf.LEFT_COUNT_KEY: 0,
                        maf.RIGHT_STATUS_KEY: "C",
                        maf.RIGHT_COUNT_KEY: 0})
  b1_tarSyr = Sequence("tarSyr1.scaffold_5923", b1_tarSyr_s,
                       8928 - 2859 - 50, 8928 - 2859, "-",
                       2859, {maf.QUALITY_META_KEY: b1_tarSyr_q,
                              maf.LEFT_STATUS_KEY: "N",
                              maf.LEFT_COUNT_KEY: 0,
                              maf.RIGHT_STATUS_KEY: "C",
                              maf.RIGHT_COUNT_KEY: 0})
  b1_mm4 = UnknownSequence("mm4.chr6", 53310102, 53310102 + 58, "+",
                           151104725 - (53310102 + 58),
                           {maf.EMPTY_ALIGNMENT_STATUS_KEY: "I"})

  # this defines a second genome alignmnet block
  b2_hg19_seq = "ccttcttttaattaattttgttaagg----gatttcctctagggccactgcacgtca"
  b2_panTro_s = "ccttcttttaattaattttgttatgg----gatttcgtctagggtcactgcacatca"
  b2_panTro_q = "99999999999999999999999999----999999099999999999999999999"
  b2_tarSyr_s = "tcttcttttaattaattttattgagggattgattccttattgggccactacacatta"
  b2_tarSyr_q = "999999899978999999999999999977989997998678865952859999899"
  b2_tupBel_s = "cct--gtttaaattactgtattg-gg----gatttcctatagggccgcttctcgtcc"
  b2_tupBel_q = "666--958759455555746366-68----656846556554745443677468565"
  b2 = "a score=31725.000000\n" +\
       "s hg19.chr22             1772 53 + 51304566 " + b2_hg19_seq + "\n" +\
       "s panTro2.chrUn          1169 53 + 58616431 " + b2_panTro_s + "\n" +\
       "q panTro2.chrUn                             " + b2_panTro_q + "\n" +\
       "i panTro2.chrUn          C 0 C 0                          " + "\n" +\
       "s tarSyr1.scaffold_5923  2909 124 -    8928 " + b2_tarSyr_s + "\n" +\
       "q tarSyr1.scaffold_5923                     " + b2_tarSyr_q + "\n" +\
       "i tarSyr1.scaffold_5923  C 0 N 0                          " + "\n" +\
       "s tupBel1.scaffold_803   33747 113 +  85889 " + b2_tupBel_s + "\n" +\
       "q tupBel1.scaffold_803                      " + b2_tupBel_q + "\n" +\
       "i tupBel1.scaffold_803 C 0 N 0              "

  # define a maf 'file' by stitching the two blocks together
  maf1 = b1 + "\n\n" + b2

  # abstract repr. of some parts of the above data.
  b2_hg19 = Sequence("hg19.chr22", b2_hg19_seq, 1772, 1825,
                     "+", 51302741)
  b2_panTro = Sequence("panTro2.chrUn", b2_panTro_s, 1169, 1169 + 53,
                       "+", 58616431 - (1169 + 53),
                       {maf.QUALITY_META_KEY: b2_panTro_q,
                        maf.LEFT_STATUS_KEY: "C",
                        maf.LEFT_COUNT_KEY: 0,
                        maf.RIGHT_STATUS_KEY: "C",
                        maf.RIGHT_COUNT_KEY: 0})
  b2_tarSyr = Sequence("tarSyr1.scaffold_5923", b2_tarSyr_s,
                       8928 - 2909 - 124, 8928 - 2909, "-",
                       2909, {maf.QUALITY_META_KEY: b2_tarSyr_q,
                              maf.LEFT_STATUS_KEY: "C",
                              maf.LEFT_COUNT_KEY: 0,
                              maf.RIGHT_STATUS_KEY: "N",
                              maf.RIGHT_COUNT_KEY: 0})


class TestGenomeAlignment(unittest.TestCase):

  """Unit tests for genome alignment IO code."""

  def setUp(self):
    """Set up a few MAF files with whole genome alignments for testing."""
    self.b1_hg19 = GATestHelper.b1_hg19
    self.b1_panTro = GATestHelper.b1_panTro
    self.b1_tarSyr = GATestHelper.b1_tarSyr
    self.b1_mm4 = GATestHelper.b1_mm4
    self.maf1 = GATestHelper.maf1
    self.b2_hg19 = GATestHelper.b2_hg19
    self.b2_panTro = GATestHelper.b2_panTro
    self.b2_tarSyr = GATestHelper.b2_tarSyr
    self.b1 = GATestHelper.b1
    self.b2 = GATestHelper.b2

  def test_genome_alignment_iterator(self):
    """Unit tests for iterating over a genome alignment in an MAF file."""
    g_iter = genome_alignment_iterator(StringIO.StringIO(self.maf1), "hg19")
    blocks = [x for x in g_iter]

    # first test that everything for regular maf stuff works properly...
    self.assertEqual(len(blocks), 2)
    self.assertEqual(blocks[0]["hg19.chr22"], self.b1_hg19)
    self.assertEqual(blocks[0]["panTro2.chrUn"], self.b1_panTro)
    self.assertEqual(blocks[0]["tarSyr1.scaffold_5923"], self.b1_tarSyr)
    self.assertEqual(blocks[0]["mm4.chr6"], self.b1_mm4)
    self.assertEqual(blocks[1]["hg19.chr22"], self.b2_hg19)
    self.assertEqual(blocks[1]["panTro2.chrUn"], self.b2_panTro)
    self.assertEqual(blocks[1]["tarSyr1.scaffold_5923"], self.b2_tarSyr)

    self.assertEqual(maf.alignment_to_maf(blocks[0]).split(), self.b1.split())
    self.assertEqual(maf.alignment_to_maf(blocks[1]).split(), self.b2.split())
    out = StringIO.StringIO()
    maf.write_maf(genome_alignment_iterator(StringIO.StringIO(self.maf1),
                                            "hg19"),
                  out)
    self.assertEqual(self.maf1.split(), out.getvalue().split())

    # now we check to see that the extra stuff granted by genome alignments
    # works...
    self.assertEqual([b.start for b in blocks], [1711, 1772])
    self.assertEqual([b.end for b in blocks], [1711 + 57, 1772 + 53])
    self.assertEqual([b.chrom for b in blocks], ["chr22", "chr22"])
    self.assertEqual([b.reference_species for b in blocks], ["hg19", "hg19"])
    self.assertEqual([b.reference_sequence_name for b in blocks],
                     ["hg19.chr22", "hg19.chr22"])

  @mock.patch('os.listdir')
  @mock.patch('os.path.isfile')
  @mock.patch('__builtin__.open')
  def test_build_genome_alignment_from_directory(self, mock_open, mock_isfile,
                                                 mock_listdir):
    """Test building a genome alignment from a whole directory of MAF files."""
    mock_listdir.return_value = ["one.maf", "two.maf", "some_sub_dir"]

    def open_side_effect(*args, **kwargs):
      if args[0] == os.path.join("the_dir", "one.maf"):
        return StringIO.StringIO(self.b1)
      elif args[0] == os.path.join("the_dir", "two.maf"):
        return StringIO.StringIO(self.b2)
      raise IOError("No such file")

    def isfile_side_effect(*args, **kwargs):
      if (args[0] == os.path.join("the_dir", "one.maf") or
          args[0] == os.path.join("the_dir", "two.maf")):
        return True
      return False

    mock_open.side_effect = open_side_effect
    mock_isfile.side_effect = isfile_side_effect

    ga = build_genome_alignment_from_directory("the_dir", "hg19")
    self.assertEqual(ga.num_blocks, 2)
    self.assertEqual(ga.get_blocks("chr22", 1711, 1720)[0]["hg19.chr22"],
                     self.b1_hg19)
    self.assertEqual(ga.get_blocks("chr22", 1770, 1780)[0]["hg19.chr22"],
                     self.b2_hg19)

  @mock.patch('__builtin__.open')
  def test_build_genome_alignment_from_file(self, mock_open):
    """Test building a genome alignment from a single MAF file."""
    def open_side_effect(*args, **kwargs):
      if args[0] == "one.maf":
        return StringIO.StringIO(self.b1 + "\n" + self.b2)
      raise IOError("No such file")

    mock_open.side_effect = open_side_effect

    ga = build_genome_alignment_from_file("one.maf", "hg19")
    self.assertEqual(ga.num_blocks, 2)
    self.assertEqual(ga.get_blocks("chr22", 1711, 1720)[0]["hg19.chr22"],
                     self.b1_hg19)
    self.assertEqual(ga.get_blocks("chr22", 1770, 1780)[0]["hg19.chr22"],
                     self.b2_hg19)

  @mock.patch('__builtin__.open')
  def test_build_genome_alignment_from_indexed_file(self, mock_open):
    """Test building GA from single MAF file which is indexed."""
    ga_in_2 = WrappedStringIO(self.maf1)
    idx_strm = StringIO.StringIO()

    # replace open with mock
    def open_side_effect(*args, **kwargs):
      if not isinstance(args[0], basestring):
        raise TypeError()
      if args[0] == "one.maf":
        return ga_in_2
      elif args[0] == "one.idx":
        return idx_strm
      raise IOError("No such file")

    mock_open.side_effect = open_side_effect

    # build an index and store it in a StringIO object
    idx_strm = _build_index(StringIO.StringIO(self.maf1), "hg19")

    ga = build_genome_alignment_from_file("one.maf", "hg19", "one.idx")
    b1_res = ga.get_blocks("chr22", 1711, 1720)[0]
    b2_res = ga.get_blocks("chr22", 1770, 1780)[0]
    num_blocks_res = ga.num_blocks
    # check that accessing num blocks and getting the blocks gave right answer,
    # but did not require a read from the MAF file.
    self.assertEqual(num_blocks_res, 2)
    self.assertEqual(b1_res.chrom, "chr22")
    self.assertEqual(b2_res.chrom, "chr22")
    self.assertEqual(b1_res.start, 1711)
    self.assertEqual(b2_res.start, 1772)
    self.assertEqual(b1_res.end, 1711 + 57)
    self.assertEqual(b2_res.end, 1772 + 53)
    self.assertEqual(ga_in_2.reads, 0)

    # check that full equality is correct, and requires loading from MAF
    self.assertEqual(b1_res["hg19.chr22"], self.b1_hg19)
    self.assertEqual(b2_res["hg19.chr22"], self.b2_hg19)
    self.assertGreater(ga_in_2.reads, 0)

  @mock.patch('os.listdir')
  @mock.patch('os.path.isfile')
  @mock.patch('__builtin__.open')
  def test_just_in_time_genome_alig(self, mock_open, mock_isfile,
                                    mock_listdir):
    """Test building a JIT genome alignment from a directory."""
    mock_listdir.return_value = ["chr22:1711-1768.maf", "chr22:1772-1825.maf",
                                 "chr22:1711-1768.idx", "chr22:1772-1825.idx",
                                 "some_sub_dir"]
    ga_in_1 = WrappedStringIO(self.b1)
    ga_in_2 = WrappedStringIO(self.b2)

    # replace open with mock
    def open_side_effect(*args, **kwargs):
      if not isinstance(args[0], basestring):
        raise TypeError()
      if args[0] == os.path.join("the_dir", "chr22:1711-1768.maf"):
        return ga_in_1
      elif args[0] == os.path.join("the_dir", "chr22:1772-1825.maf"):
        return ga_in_2
      elif args[0] == os.path.join("the_dir", "chr22:1711-1768.idx"):
        x = _build_index(StringIO.StringIO(self.b1), "hg19")
        x.seek(0)
        return x
      elif args[0] == os.path.join("the_dir", "chr22:1772-1825.idx"):
        x = _build_index(StringIO.StringIO(self.b2), "hg19")
        x.seek(0)
        return x
      raise IOError("No such file: " + args[0])

    def isfile_side_effect(*args, **kwargs):
      fn = args[0].strip()
      if (fn == os.path.join("the_dir", "chr22:1711-1768.maf") or
          fn == os.path.join("the_dir", "chr22:1772-1825.maf") or
          fn == os.path.join("the_dir", "chr22:1772-1825.idx") or
          fn == os.path.join("the_dir", "chr22:1711-1768.idx")):
        return True
      if fn == os.path.join("the_dir", "some_sub_dir"):
        return False
      raise IOError("No such file: " + fn)

    mock_open.side_effect = open_side_effect
    mock_isfile.side_effect = isfile_side_effect

    # b1 -> chr22:1711-1768
    # b2 -> chr22:1772-1825
    ga = load_just_in_time_genome_alignment("the_dir", "hg19",
                                            extensions=["maf"],
                                            index_exts=["idx"])
    b1_res = ga.get_blocks("chr22", 1711, 1720)[0]
    b2_res = ga.get_blocks("chr22", 1770, 1780)[0]

    # check that accessing num blocks and getting the blocks gave right answer,
    # but did not require a read from the MAF file.
    self.assertEqual(b1_res.chrom, "chr22")
    self.assertEqual(b2_res.chrom, "chr22")
    self.assertEqual(b1_res.start, 1711)
    self.assertEqual(b2_res.start, 1772)
    self.assertEqual(b1_res.end, 1711 + 57)
    self.assertEqual(b2_res.end, 1772 + 53)
    self.assertEqual(ga_in_1.reads, 0)
    self.assertEqual(ga_in_2.reads, 0)

    # check that full equality is correct, and requires loading from MAF
    self.assertEqual(b1_res["hg19.chr22"], self.b1_hg19)
    self.assertEqual(b2_res["hg19.chr22"], self.b2_hg19)
    self.assertGreater(ga_in_1.reads, 0)
    self.assertGreater(ga_in_2.reads, 0)


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
