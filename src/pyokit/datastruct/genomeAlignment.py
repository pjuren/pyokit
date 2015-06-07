"""
  Date of Creation: 26th May 2015

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


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class GenomeAlignmentError(Exception):
  """
  Class representing errors that occur when manipulating whole-genome
  alignments.
  """
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


###############################################################################
#                              HELPER FUNCTIONS                               #
###############################################################################

def _build_trees_by_chrom(blocks):
  """
  Construct a set of interval trees from an iterable of genome alignment blocks
  such that each tree represents a chromosome.
  :return: a dictionary indexed by chromosome name where each entry is an
           interval tree
  """
  by_chrom = {}
  for b in blocks:
    if b.chrom not in by_chrom:
      by_chrom[b.chrom] = []
    by_chrom[b.chrom].append(b)
  res = {}
  for c in by_chrom:
    res[c] = IntervalTree(by_chrom[c])
  return res


###############################################################################
#                        GENOME ALIGNMENT BLOCK CLASS                         #
###############################################################################

class GenomeAlignmentBlock(MultipleSequenceAlignment):
  def __init__(self, sequences, reference_species, meta_data=None):
    """
    :param sequences:     the sequences that align to this block of the
                          reference. Note that the reference must be named
                          assembly.chrom, e.g. hg19.chr13.
    :param reference_seq: the name of the reference species assembly,
                          e.g. hg19. Note that this must match up with the name
                          given to the sequence in the sequences list.
    """
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
    return self[self.reference_sequence_name].start

  @property
  def end(self):
    return self[self.reference_sequence_name].end


###############################################################################
#                   ON-DEMAND GENOME ALIGNMENT BLOCK CLASS                    #
###############################################################################

@decorate_all_methods_and_properties(just_in_time_method,
                                     just_in_time_property)
class JustInTimeGenomeAlignmentBlock(GenomeAlignmentBlock):
  """
  A genome alignment block that is loaded just-in-time from some factory
  object; a common pattern would be to provide an IndexedFile as the factory
  object

  :param factory: any object that implements the subscript operator such that
                  it accept the key as a unique identifier and returns the
                  alignment block
  :param key:     any object which uniquely identifies this pariwise alignment
                  to the factory; i.e. a hash key
  """
  def __init__(self, factory, key):
    self.factory = factory
    self.key = key
    self.item = None


###############################################################################
#                           GNEOME ALIGNMENT CLASS                            #
###############################################################################

class GenomeAlignment(object):
  """
  :param blocks: ???
  """

  def __init__(self, blocks):
    """
    Constructor for GenomeAlignment objects; see class docstring for
    description of parameters.
    """
    # a genome alignment stores blocks internally using one an interval tree
    # for each chromosome
    self.block_trees = _build_trees_by_chrom(blocks)
    self.num_blocks = len(blocks)

  def get_blocks(self, chrom, start, end):
    """
    :return: the alignment blocks that overlap a given genomic interval;
             potentially none, in which case the empty list is returned.
    """
    if chrom not in self.block_trees:
      return []
    return self.block_trees[chrom].intersectingInterval(start, end)


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################
class TestGenomeAlignment(unittest.TestCase):
  pass


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
