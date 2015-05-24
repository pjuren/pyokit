"""
  Date of Creation: 22nd May 2015

  Description:   Classes and for representing pairwise and multiple sequence
                 alignments.

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
import unittest
import StringIO

# pyokit imports
from pyokit.datastruct.sequence import Sequence
from pyokit.datastruct.multipleAlignment import MultipleSequenceAlignment


###############################################################################
#                                 CONSTANTS                                   #
###############################################################################

# the following are the characters that mark the start of certain line types
# in MAF files. See the header of the MAF iterator below for more details
A_LINE = "a"
S_LINE = "s"
I_LINE = "i"
E_LINE = "e"
Q_LINE = "q"


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class MAFError(Exception):
  """
  Class representing errors that occur when manipulating pairwise or multiple
  alignment objects
  """
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


###############################################################################
#                                 ITERATORS                                   #
###############################################################################

def maf_iterator(fn):
  """
  MAF files are arranged in blocks. Each block is a multiple alignment. Within
  a block, the first character of a line indicates what kind of line it is:

  a -- gives key-value pair meta data for the block (should only be one of
       these lines, and it should be the first line in the block)
  s -- a sequence line; has 6 fields in addition to 's':
          * source sequence,
          * start coordinate of the sequence, Zero-based. If strand is -'ve,
            this is relative to the start of the reverse complement of the
            source sequence.
          * ungapped length of the sequence
          * strand
          * src size -- the full length of the source sequence
          * the sequence itself
  i -- always come after s lines, and contain information about the context of
       the sequence. Five fields are given, not counting the 'i'
          * source sequence (must match s line before this)
          * left status (see below)
          * left count; num of bases in source sequence between start of the
            block and end of previous block (0 if this is the first)
          * right status (see below)
          * left count; num of bases in source after end of this block before
            start of next
       status (left/right) is a single char and can be:
          * C -- the sequence before or after is contiguous with this block.
          * I -- there are bases between the bases in this block and the one
                 before or after it.
          * N -- this is the first sequence from this src chrom or scaffold.
          * n -- this is the first sequence from this src chrom or scaffold but
                 it is bridged by another alignment from a different chrom or
                 scaffold.
          * M -- there is missing data before or after this block (Ns in the
                 sequence).
          * T -- the sequence in this block has been used before in a previous
                 block (likely a tandem duplication)
  e -- indicates that there is no aligning sequence for a species, but that
       there are blocks before and after this one that do align. Format is the
       same as the 's' lines, but the start and length indicate the start and
       size of the non-aligning region in the sequence and the sequence is
       replaced with a status character:
       C -- the sequence before and after is contiguous implying that this
            region was either deleted in the source or inserted in the
            reference sequence.
       I -- there are non-aligning bases in the source species between
            chained alignment blocks before and after this block.
       M -- there are non-aligning bases in the source and more than 90% of
            them are Ns in the source. The browser shows a pale yellow bar.
       n -- there are non-aligning bases in the source and the next aligning
            block starts in a new chromosome or scaffold that is bridged
            by a chain between still other blocks.
  q -- quality information about an aligned base in a species. Two fields after
       the 'q': the source name and a single digit for each nucleotide in its
       sequence (0-9 or F, or - to indicate a gap).
  """

  try:
    fh = open(fn)
  except (TypeError):
    fh = fn

  sequences = []
  meta_data = {}
  for line in fh:
    line = line.strip()
    if line == "":
      continue
    parts = line.split()
    line_type = parts[0].strip()
    if line_type == A_LINE:
      if sequences != []:
        yield MultipleSequenceAlignment(sequences, meta_data)
      sequences = []
      meta_data = {}
      for i in range(1, len(parts)):
        assert(parts[i].count("=") == 1)
        piv = parts[i].find("=")
        meta_data[parts[i][:piv]] = parts[i][piv + 1:]
    elif line_type == S_LINE:
      strand = parts[4]
      total_seq_len = int(parts[5])
      start = int(parts[2]) if strand == "+" else total_seq_len - int(parts[2])
      end = start + int(parts[3])
      remain = total_seq_len - end
      sequences.append(Sequence(parts[1], parts[6], start, end,
                                strand, remain))
    elif line_type == I_LINE:
      pass
    elif line_type == E_LINE:
      pass
    elif line_type == Q_LINE:
      pass
    else:
      raise MAFError("Unknown type of MAF line: " + line)

  # don't forget to yield the final block
  if sequences != []:
    yield MultipleSequenceAlignment(sequences, meta_data)


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################

class TestMAF(unittest.TestCase):
  def setUp(self):
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
         "i tupBel1.scaffold_803   I 1 C 0            "

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
    self.maf1 = b1 + "\n\n" + b2

  def test_maf_iterator(self):
    in_file_contents = StringIO.StringIO(self.maf1)
    blocks = [x for x in maf_iterator(in_file_contents)]
    self.assertEqual(len(blocks), 2)

###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
