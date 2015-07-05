"""
Date of Creation: 26th May 2015.

Description:   Whole-genome alignment a reference genome and one or more
               others

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

# standard imports
import unittest

# pyokit imports
from pyokit.datastruct.multipleAlignment import MultipleSequenceAlignment
from pyokit.datastruct.multipleAlignment import MissingSequenceHandler
from pyokit.datastruct.intervalTree import IntervalTree
from pyokit.util.meta import decorate_all_methods
from pyokit.util.meta import just_in_time_method
from pyokit.datastruct.sequence import Sequence
from pyokit.datastruct.genomicInterval import GenomicInterval


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class GenomeAlignmentError(Exception):

  """Errors that occur when manipulating whole-genome alignments."""

  def __init__(self, msg):
    """
    Constructor for GenomeAlignmentErrors.

    :param msg: message to display
    """
    self.value = msg

  def __str__(self):
    """:return: string representation of error."""
    return repr(self.value)


class NoSuchAlignmentColumnError(GenomeAlignmentError):

  """Raised when a column is requested that doesn't exist in the alignment."""

  def __init__(self, msg):
    """
    Constructor for GenomeAlignmentErrors.

    :param msg: message to display
    """
    self.value = msg


class NoUniqueColumnError(GenomeAlignmentError):

  """Raised when no unique alignment present (i.e. ambiguous)."""

  def __init__(self, msg):
    """
    Constructor for GenomeAlignmentErrors.

    :param msg: message to display
    """
    self.value = msg


###############################################################################
#                              HELPER FUNCTIONS                               #
###############################################################################

def _build_trees_by_chrom(blocks):
  """
  Construct set of interval trees from an iterable of genome alignment blocks.

  :return: a dictionary indexed by chromosome name where each entry is an
           interval tree for that chromosome.
  """
  by_chrom = {}
  for b in blocks:
    if b.chrom not in by_chrom:
      by_chrom[b.chrom] = []
    by_chrom[b.chrom].append(b)
  res = {}
  for c in by_chrom:
    res[c] = IntervalTree(by_chrom[c], openEnded=True)
  return res


###############################################################################
#                        GENOME ALIGNMENT BLOCK CLASS                         #
###############################################################################

class GenomeAlignmentBlock(MultipleSequenceAlignment):

  """
  Class representing a single block in a whole-genome alignment.

  Blocks are segments of the reference genome that are aligned to one or more
  or the query genomes.

  :param sequences:     the sequences that align to this block of the
                        reference. Note that the reference must be named
                        assembly.chrom, e.g. hg19.chr13.
  :param reference_seq: the name of the reference species assembly,
                        e.g. hg19. Note that this must match up with the name
                        given to the sequence in the sequences list.
  """

  def __init__(self, sequences, reference_species, meta_data=None):
    """Constructor; see class docstring for parameter details."""
    MultipleSequenceAlignment.__init__(self, sequences, meta_data)

    ref_index = None
    for i in range(0, len(sequences)):
      if sequences[i].name[0:len(reference_species)] == reference_species:
        if ref_index is not None:
          raise GenomeAlignmentError("multiple sequences in block match " +
                                     "reference species name: " +
                                     sequences[i] + " and " +
                                     sequences[ref_index])
        ref_index = i
    self.reference_sequence_name = sequences[ref_index].name
    name_parts = self.reference_sequence_name.split(".")
    if len(name_parts) != 2:
      raise GenomeAlignmentError("invalid name for genome alignment block: " +
                                 sequences[ref_index].name + "; must be " +
                                 "formatted assembly.chrom, e.g. hg19.chr1")
    self._chrom = name_parts[1]
    self.reference_species = reference_species

  @property
  def chrom(self):
    """:return the chromosome the block is on."""
    return self._chrom

  @property
  def start(self):
    """:return: the start of the block on the reference genome sequence."""
    return self[self.reference_sequence_name].start

  @property
  def end(self):
    """:return: the end of the block on the reference genome sequence."""
    return self[self.reference_sequence_name].end

  def get_column_absolute(self, position):
    """
    return a column from the block as dictionary indexed by seq. name.

    :param position: the index to extract from the block; must be absolute
                     coordinates (i.e. between self.start and self.end, not
                     inclusive of the end).
    :return: dictionary where keys are sequence names and values are
             nucleotides (raw strings).
    """
    if position < self.start or position >= self.end:
      raise ValueError("getting column at genomic locus " + self._chrom + " " +
                       str(position) + " failed; locus is outside of genome " +
                       "alignment block")

    rel_coord = self.sequence_to_alignment_coords(self.reference_sequence_name,
                                                  position, position + 1)
    assert(len(rel_coord) == 1)
    rel_start, rel_end = rel_coord[0]
    assert(rel_end == rel_start + 1)
    return self.get_column(rel_start, MissingSequenceHandler.TREAT_AS_ALL_GAPS)


###############################################################################
#                   ON-DEMAND GENOME ALIGNMENT BLOCK CLASS                    #
###############################################################################

@decorate_all_methods(just_in_time_method)
class JustInTimeGenomeAlignmentBlock(GenomeAlignmentBlock):

  """
  Genome alignment block loaded just-in-time from some factory object.

  A common pattern would be to provide an IndexedFile as the factory object.
  The start, end and chrom are provided in the key, which should be
  constructed using the static method build_hash. Accessing these will not
  trigger loading the whole algnment, but calls to any other methods or
  properties will cause a full load.

  :param factory: any object that implements the subscript operator such that
                  it accept the key as a unique identifier and returns the
                  alignment block. Probably an IndexedFile, but doesn't have
                  to be.
  :param key:     a hash for this genome alignment; build it using the static
                  method build_hash.
  """

  HASH_SEP = "\t"

  def __init__(self, factory, key):
    """Constructor; see class docstring for parameter details."""
    self.factory = factory
    self.key = key
    self.item = None

    self._chrom, self._start, self._end =\
        JustInTimeGenomeAlignmentBlock.parse_hash(key)

  @property
  def chrom(self):
    """:return: the chromosome this block is one in the reference seq."""
    if self.item is None:
      return self._chrom
    else:
      if self._chrom != self.item.chrom:
        raise IndexError("Chromosome mismatch between index and genome " +
                         "alignment block")
      return self.item.chrom

  @property
  def start(self):
    """:return: the genomic start of this alignment block in ref seq."""
    if self.item is None:
      return self._start
    else:
      if self._start != self.item.start:
        raise IndexError("Start mismatch between index and genome " +
                         "alignment block")
      return self.item.start

  @property
  def end(self):
    """:return: the genomic end of this alignment block in ref seq."""
    if self.item is None:
      return self._end
    else:
      if self._end != self.item.end:
        raise IndexError("End mismatch between index and genome " +
                         "alignment block")
      return self.item.end

  @staticmethod
  def build_hash(alig):
    """Build an index hash for an alignment at a given genomic locus."""
    return JustInTimeGenomeAlignmentBlock.build_hash_raw(alig.chrom,
                                                         alig.start, alig.end)

  @staticmethod
  def build_hash_raw(chrom, start, end):
    """Build a hash string for an alignment from raw chrom, start, end."""
    return chrom + "\t" + str(start) + "\t" + str(end)

  @staticmethod
  def parse_hash(h):
    """Extract the chrom, start and end coordinates from a hash value."""
    chrom, start, end = h.split(JustInTimeGenomeAlignmentBlock.HASH_SEP)
    start = int(start)
    end = int(end)
    return chrom, start, end


###############################################################################
#                           GNEOME ALIGNMENT CLASS                            #
###############################################################################

class GenomeAlignment(object):

  """
  A genome alignment is a aggregation of genome alignment blocks.

  :param blocks: genome alignment blocks to build the alignment from.
  """

  def __init__(self, blocks):
    """Constructor; see class docstring for description of parameters."""
    # a genome alignment stores blocks internally using one an interval tree
    # for each chromosome
    self.block_trees = _build_trees_by_chrom(blocks)
    self.num_blocks = len(blocks)

  def get_blocks(self, chrom, start, end):
    """
    Get any blocks in this alignment that overlap the given location.

    :return: the alignment blocks that overlap a given genomic interval;
             potentially none, in which case the empty list is returned.
    """
    if chrom not in self.block_trees:
      return []
    return self.block_trees[chrom].intersectingInterval(start, end)

  def get_column(self, chrom, position):
    """Get the alignment column at the specified chromosome and position."""
    blocks = self.get_blocks(chrom, position, position + 1)
    if len(blocks) == 0:
      raise NoSuchAlignmentColumnError("Request for column on chrom " +
                                       chrom + " at position " +
                                       str(position) + " not possible; " +
                                       "genome alignment not defined at " +
                                       "that locus.")
    if len(blocks) > 1:
      raise NoUniqueColumnError("Request for column on chrom " + chrom +
                                " at position " + str(position) + "not " +
                                "possible; ambiguous alignment of that locus.")

    return blocks[0].get_column_absolute(position)


###############################################################################
#                    JUST-IN-TIME GNEOME ALIGNMENT CLASS                      #
###############################################################################

class JITGenomeAlignmentKeyInterval(object):
  def __init__(self, interval, key):
    self.interval = interval
    self.key = key


class JustInTimeGenomeAlignment(GenomeAlignment):

  """
  A genome alignment stored on-dsik over many alignment and index files.

  Only one file will be loaded at a time, and loading will be done
  just-in-time to satisfy lookups

  :param whole_chrom_files:   a dictionary indexed by chrom name where each
                              value is a hash key that uniquely identifies the
                              corresponding genome alignment to the factory
  :param partial_chrom_files: a dictionary indexed by tuple of (chrom, start,
                              end). each value is a hash key that uniquely
                              identifies the corresponding genome alignment to
                              the factory.
  :param: factory             A callable that accepts the hash keys from
                              whole_chrom_files and partial_chrom_files and
                              returns the corresponding genome alignment.
  """

  def __init__(self, whole_chrom_files, partial_chrom_files, factory):
    """Constructor; see class docsstring for param details."""
    self.current = None
    self.whole_chrom_files = whole_chrom_files
    self.partial_trees = {}
    by_chrom = {}
    for chrom, start, end in partial_chrom_files:
      k = (chrom, start, end)
      v = partial_chrom_files[k]
      if chrom in whole_chrom_files:
        raise GenomeAlignmentError("Oops")
      if chrom not in by_chrom:
        by_chrom[chrom] = []
      interval = GenomicInterval(chrom, start, end)
      by_chrom[chrom].append(JITGenomeAlignmentKeyInterval(interval, v))
    for chrom in by_chrom:
      self.partial_trees[chrom] = IntervalTree(by_chrom[chrom])
    for chrom, start, end in partial_chrom_files:
      hits = self.partial_trees[chrom].intersecting_interval(start, end)
      if len(hits) != 1:
        raise GenomeAlignmentError("Oops")

  def __get_keys(self, chrom, start, end=None):
    """
    Get the hash keys in whole/partial_chrom_files that overlap the interval.

    :return: list of hash keys.
    """
    keys = []
    if chrom in self.whole_chrom_files:
      keys.append(self.whole_chrom_files[chrom])
    if chrom in self.partial_trees:
      if end is not None:
        hits = self.partial_trees[chrom].intersectingInterval(start, end)
      else:
        hits = self.partial_trees[chrom].intersectingPoint(start, end)
      for hit in hits:
        keys.append(hit.key)
    return keys

  def __switch_alig(self, key):
    self.current = self.factory[key]

  def get_blocks(self, chrom, start, end):
    """
    Get any blocks in this alignment that overlap the given location.

    :return: the alignment blocks that overlap a given genomic interval;
             potentially none, in which case the empty list is returned.
    """
    blocks = []
    for k in self.__get_keys(chrom, start, end):
      self.__switch_alig(k)
      blocks += self.current.get_blocks(chrom, start, end)
    return blocks


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################
class TestGenomeAlignmentDS(unittest.TestCase):

  """Unit tests for GenomeAlignment module."""

  def setUp(self):
    """Set up some genome alignment blocks and genome alignments for tests."""
    self.block1 = GenomeAlignmentBlock([Sequence("s1.c1", "TCTCGC-A", 11, 18),
                                        Sequence("s2.c1", "ACTGGC--", 25, 31),
                                        Sequence("s3.c2", "ACTGCCTA", 5, 13),
                                        Sequence("s4.c1", "ACT-GCTA", 58, 65)],
                                       "s1")
    self.block2 = GenomeAlignmentBlock([Sequence("s1.c2", "C------G", 21, 23),
                                        Sequence("s2.c2", "CGGTCAGG", 85, 94),
                                        Sequence("s3.c2", "-GGTC-GG", 1, 7),
                                        Sequence("s4.c3", "-GGCCAGG", 3, 11)],
                                       "s1")
    # this block defines an ambiguous alignment of part of block1
    self.block3 = GenomeAlignmentBlock([Sequence("s1.c1", "GCACGCT", 15, 22),
                                        Sequence("s2.c8", "GCAC-CT", 25, 31),
                                        Sequence("s3.c8", "GC-CGCT", 5, 13),
                                        Sequence("s4.c8", "GC--GCT", 58, 65)],
                                       "s1")
    self.ga1 = GenomeAlignment([self.block1, self.block2])
    self.ga2 = GenomeAlignment([self.block1, self.block2, self.block3])

  def test_block_get_column(self):
    """Test getting a single column from a block."""
    exp1 = {"s1.c1": "T", "s2.c1": "A", "s3.c2": "A", "s4.c1": "A"}
    self.assertEqual(self.block1.get_column_absolute(11), exp1)
    exp2 = {"s1.c1": "A", "s2.c1": "-", "s3.c2": "A", "s4.c1": "A"}
    self.assertEqual(self.block1.get_column_absolute(17), exp2)
    self.assertRaises(ValueError, self.block1.get_column_absolute, 18)
    self.assertRaises(ValueError, self.block1.get_column_absolute, 10)
    exp3 = {"s1.c2": "G", "s2.c2": "G", "s3.c2": "G", "s4.c3": "G"}
    self.assertEqual(self.block2.get_column_absolute(22), exp3)

  def test_genome_alig_get_col(self):
    """Test getting a single column from a genome alignment."""
    exp1 = {"s1.c1": "T", "s2.c1": "A", "s3.c2": "A", "s4.c1": "A"}
    self.assertEqual(self.ga1.get_column("c1", 11), exp1)
    exp2 = {"s1.c2": "C", "s2.c2": "C", "s3.c2": "-", "s4.c3": "-"}
    self.assertEqual(self.ga1.get_column("c2", 21), exp2)
    self.assertRaises(NoUniqueColumnError, self.ga2.get_column, "c1", 15)
    self.assertRaises(NoSuchAlignmentColumnError, self.ga1.get_column,
                      "c1", 10)
    self.assertRaises(NoSuchAlignmentColumnError, self.ga1.get_column,
                      "c3", 15)

  def test_jit_genome_alig_block(self):
    """Test wrapping a block with JIT wrapper."""
    b1_hash = JustInTimeGenomeAlignmentBlock.build_hash_raw("chr1", 11, 18)
    factory = {b1_hash: self.block1}
    b1_jit = JustInTimeGenomeAlignmentBlock(factory, b1_hash)
    expect = {"s1.c1": "T", "s2.c1": "A", "s3.c2": "A", "s4.c1": "A"}
    self.assertEqual(b1_jit.get_column_absolute(11), expect)

  def test_jit_fail_on_mismatch(self):
    """Test that JIT genome alig blocks fail when index doesn't match."""
    def chrom_prop_wrapper(b):
      return b.chrom

    def start_prop_wrapper(b):
      return b.start

    def end_prop_wrapper(b):
      return b.end

    # for our sanity, check that it works when everything is given correctly.
    b1_hash = JustInTimeGenomeAlignmentBlock.build_hash_raw("c1", 11, 18)
    factory = {b1_hash: self.block1}
    b1_jit = JustInTimeGenomeAlignmentBlock(factory, b1_hash)
    self.assertEqual(b1_jit.chrom, "c1")
    self.assertEqual(b1_jit.start, 11)
    self.assertEqual(b1_jit.end, 18)
    b1_jit.get_column_absolute(11)  # this triggers a full load from factory
    self.assertEqual(b1_jit.chrom, "c1")
    self.assertEqual(b1_jit.start, 11)
    self.assertEqual(b1_jit.end, 18)

    # wrong chrom
    b1_hash = JustInTimeGenomeAlignmentBlock.build_hash_raw("chr2", 11, 18)
    factory = {b1_hash: self.block1}
    b1_jit = JustInTimeGenomeAlignmentBlock(factory, b1_hash)
    self.assertEqual(b1_jit.chrom, "chr2")
    b1_jit.get_column_absolute(11)  # this triggers a full load from factory
    self.assertRaises(IndexError, chrom_prop_wrapper, b1_jit)
    self.assertEqual(b1_jit.start, 11)  # should still be okay though
    self.assertEqual(b1_jit.end, 18)    # should still be okay though

    # wrong start
    b1_hash = JustInTimeGenomeAlignmentBlock.build_hash_raw("c1", 10, 18)
    factory = {b1_hash: self.block1}
    b1_jit = JustInTimeGenomeAlignmentBlock(factory, b1_hash)
    self.assertEqual(b1_jit.chrom, "c1")
    b1_jit.get_column_absolute(11)  # this triggers a full load from factory
    self.assertRaises(IndexError, start_prop_wrapper, b1_jit)
    self.assertEqual(b1_jit.chrom, "c1")  # should still be okay though
    self.assertEqual(b1_jit.end, 18)    # should still be okay though

    # wrong end
    b1_hash = JustInTimeGenomeAlignmentBlock.build_hash_raw("c1", 11, 17)
    factory = {b1_hash: self.block1}
    b1_jit = JustInTimeGenomeAlignmentBlock(factory, b1_hash)
    self.assertEqual(b1_jit.chrom, "c1")
    b1_jit.get_column_absolute(11)  # this triggers a full load from factory
    self.assertRaises(IndexError, end_prop_wrapper, b1_jit)
    self.assertEqual(b1_jit.start, 11)  # should still be okay though
    self.assertEqual(b1_jit.chrom, "c1")    # should still be okay though

  def test_jit_genome_alig_block_no_unneeded_loads(self):
    """Test using a JIT genome alig block without causing a full load."""
    class SimpleFactory:
      def __getitem__(self, k):
        raise ValueError()

    b1_hash = JustInTimeGenomeAlignmentBlock.build_hash_raw("chr1", 11, 18)
    b1_jit = JustInTimeGenomeAlignmentBlock(SimpleFactory(), b1_hash)

    # these should be fine
    self.assertEqual(b1_jit.chrom, "chr1")
    self.assertEqual(b1_jit.start, 11)
    self.assertEqual(b1_jit.end, 18)

    # these should fail
    self.assertRaises(ValueError, b1_jit.get_column_absolute, 11)

###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
