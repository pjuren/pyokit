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
from pyokit.datastruct.intervalTree import IntervalTree
from pyokit.util.meta import decorate_all_methods_and_properties
from pyokit.util.meta import just_in_time_method
from pyokit.util.meta import just_in_time_property
from pyokit.datastruct.sequence import Sequence


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
    self.chrom = name_parts[1]
    self.reference_species = reference_species

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
      raise ValueError("getting column at genomic locus " + self.chrom + " " +
                       str(position) + " failed; locus is outside of genome " +
                       "alignment block")

    rel_coord = self.sequence_to_alignment_coords(self.reference_sequence_name,
                                                  position, position + 1)
    assert(len(rel_coord) == 1)
    rel_start, rel_end = rel_coord[0]
    assert(rel_end == rel_start + 1)
    return self.get_column(rel_start)


###############################################################################
#                   ON-DEMAND GENOME ALIGNMENT BLOCK CLASS                    #
###############################################################################

@decorate_all_methods_and_properties(just_in_time_method,
                                     just_in_time_property)
class JustInTimeGenomeAlignmentBlock(GenomeAlignmentBlock):

  """
  Genome alignment block loaded just-in-time from some factory object.

  A common pattern would be to provide an IndexedFile as the factory object.
  Any access to the genome alignment block's methods will trigger loading
  the full alignment from the factory, otherwise it is never loaded.

  :param factory: any object that implements the subscript operator such that
                  it accept the key as a unique identifier and returns the
                  alignment block
  :param key:     any object which uniquely identifies this pariwise alignment
                  to the factory; i.e. a hash key
  """

  def __init__(self, factory, key):
    """Constructor; see class docstring for parameter details."""
    self.factory = factory
    self.key = key
    self.item = None


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

###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
