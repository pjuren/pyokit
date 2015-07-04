"""
Date of Creation: 22nd May 2015.

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
from pyokit.io.ioError import PyokitIOError
from pyokit.datastruct.sequence import Sequence
from pyokit.datastruct.sequence import UnknownSequence
from pyokit.datastruct.multipleAlignment import MultipleSequenceAlignment


###############################################################################
#                                 CONSTANTS                                   #
###############################################################################
# special meta-data key for alignments that stores the order the sequences
# were presented in, so the output can reproduce it.
SEQ_ORDER_KEY = "SEQ_ORDER_KEY"

# the following are the characters that mark the start of certain line types
# in MAF files. See the header of the MAF iterator below for more details
A_LINE = "a"
S_LINE = "s"
I_LINE = "i"
E_LINE = "e"
Q_LINE = "q"

# keys for placing meta data into Sequence meta-data dictionaries
QUALITY_META_KEY = "QL"
EMPTY_ALIGNMENT_STATUS_KEY = "EMPTY_STATUS"
LEFT_STATUS_KEY = "LEFT_STATUS"
LEFT_COUNT_KEY = "LEFT_COUNT"
RIGHT_STATUS_KEY = "RIGHT_STATUS"
RIGHT_COUNT_KEY = "RIGHT_COUNT"


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class MAFError(PyokitIOError):

  """Errors that occur when performing IO on .maf format files."""

  def __init__(self, msg):
    """Constructor for MAFErrors."""
    self.value = msg

  def __str__(self):
    """:return: string representation of this MAFError."""
    return repr(self.value)


###############################################################################
#                             HELPER FUNCTIONS                                #
###############################################################################

def merge_dictionaries(a, b):
  """Merge two dictionaries; duplicate keys get value from b."""
  res = {}
  for k in a:
    res[k] = a[k]
  for k in b:
    res[k] = b[k]
  return res


def __build_sequence(parts):
  """Build a sequence object using the pre-tokenized parts from a MAF line.

  s -- a sequence line; has 6 fields in addition to 's':
          * source sequence,
          * start coord. of seq., zero-based. If -'ve strand, rel to start of
            rev. comp.
          * ungapped length of the sequence
          * strand
          * src size -- the full length of the source sequence
          * the sequence itself
  """
  strand = parts[4]
  seq_length = int(parts[3])
  total_seq_len = int(parts[5])
  start = (int(parts[2]) if strand == "+"
           else total_seq_len - int(parts[2]) - seq_length)
  end = start + seq_length
  remain = total_seq_len - end
  return Sequence(parts[1], parts[6], start, end, strand, remain)


def __build_unknown_sequence(parts):
  """Build unknown seq (e lines) using pre-tokenized parts from MAF line.

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
  """
  strand = parts[4]
  seq_length = int(parts[3])
  total_seq_len = int(parts[5])
  start = (int(parts[2]) if strand == "+"
           else total_seq_len - int(parts[2]) - seq_length)
  end = start + seq_length
  remain = total_seq_len - end
  return UnknownSequence(parts[1], start, end, strand, remain,
                         {EMPTY_ALIGNMENT_STATUS_KEY: parts[6]})


def __annotate_sequence_with_context(seq, i_line_parts):
  """Extract meta data from pre-tokenized maf i-line and populate sequence.

  i -- always come after s lines, and contain information about the context of
       the sequence. Five fields are given, not counting the 'i'
          * source sequence (must match s line before this)
          * left status (see below)
          * left count; num of bases in source sequence between start of the
            block and end of previous block (0 if this is the first)
          * right status (see below)
          * right count; num of bases in source after end of this block before
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
  """
  if i_line_parts[1] != seq.name:
    raise MAFError("Trying to populate meta data for sequence " + seq.name +
                   " with i-line information for " +
                   str(i_line_parts[1]) + "; maflormed MAF file?")
  if len(i_line_parts) != 6:
    raise MAFError("i-line with " + str(len(i_line_parts)) + "; expected 6.")
  seq.meta_data[LEFT_STATUS_KEY] = i_line_parts[2]
  seq.meta_data[LEFT_COUNT_KEY] = int(i_line_parts[3])
  seq.meta_data[RIGHT_STATUS_KEY] = i_line_parts[4]
  seq.meta_data[RIGHT_COUNT_KEY] = int(i_line_parts[5])


def __annotate_sequence_with_quality(seq, q_line_parts):
  """Extract meta data from pre-tokenized maf q-line and populate sequence.

  q -- quality information about an aligned base in a species. Two fields after
       the 'q': the source name and a single digit for each nucleotide in its
       sequence (0-9 or F, or - to indicate a gap).
  """
  if q_line_parts[1] != seq.name:
    raise MAFError("trying to populate meta data for sequence " + seq.name +
                   " with q-line information for " +
                   str(q_line_parts[1]) + "; maflormed MAF file?")
  if len(q_line_parts[2]) != len(seq):
    raise MAFError("trying to populate quality meta data for sequence with " +
                   "length " + str(len(seq)) + " using quality line with " +
                   "length " + str(len(q_line_parts[2])) + "; malformed " +
                   "MAF file?")
  seq.meta_data[QUALITY_META_KEY] = q_line_parts[2]


###############################################################################
#                                 ITERATORS                                   #
###############################################################################

def maf_iterator(fn, index_friendly=False,
                 yield_class=MultipleSequenceAlignment, yield_kw_args={},
                 verbose=False):
  """
  Iterate of MAF format file and yield <yield_class> objects for each block.

  MAF files are arranged in blocks. Each block is a multiple alignment. Within
  a block, the first character of a line indicates what kind of line it is:

  a -- key-value pair meta data for block; one per block, should be first line
  s -- a sequence line; see __build_sequence() for details.
  i -- always come after s lines, and contain information about the context of
       the sequence; see __annotate_sequence_with_context() for details.
  e -- indicates that there is no aligning sequence for a species, but that
       there are blocks before and after this one that do align. See
       __build_unknown_sequence() for details.
  q -- quality information about an aligned base in a species. See
       __annotate_sequence_with_quality()

  :param yield_class:   yield objects returned by this function/constructor;
                        must accept key-word args of 'sequences' and
                        'meta_data' which are list and dictionary respectively.
                        Default is MultipleSequenceAlignment
  :param yield_kw_args: extra keyword args to pass to 'yield_class'
  """

  try:
    fh = open(fn)
  except (TypeError):
    fh = fn
  if index_friendly:
    fh = iter(fh.readline, '')

  sequences = []
  meta_data = {}
  for line in fh:
    line = line.strip()
    if line == "" or line[0] == "#":
      continue
    parts = line.split()
    line_type = parts[0].strip()
    if line_type == A_LINE:
      if sequences != []:
        meta_data[SEQ_ORDER_KEY] = [s.name for s in sequences]
        kw_args = merge_dictionaries({"sequences": sequences,
                                      "meta_data": meta_data}, yield_kw_args)
        yield yield_class(**kw_args)
      sequences = []
      meta_data = {}
      for i in range(1, len(parts)):
        assert(parts[i].count("=") == 1)
        piv = parts[i].find("=")
        meta_data[parts[i][:piv]] = parts[i][piv + 1:]
    elif line_type == S_LINE:
      sequences.append(__build_sequence(parts))
    elif line_type == I_LINE:
      if len(sequences) < 1:
        raise MAFError("found information line with no preceeding sequence " +
                       "in block")
      __annotate_sequence_with_context(sequences[-1], parts)
    elif line_type == E_LINE:
      sequences.append(__build_unknown_sequence(parts))
    elif line_type == Q_LINE:
      if len(sequences) < 1:
        raise MAFError("found quality line with no preceeding sequence in " +
                       "block")
      __annotate_sequence_with_quality(sequences[-1], parts)
    else:
      raise MAFError("Unknown type of MAF line: " + line)

  # don't forget to yield the final block
  if sequences != []:
    meta_data[SEQ_ORDER_KEY] = [s.name for s in sequences]
    kw_args = merge_dictionaries({"sequences": sequences,
                                  "meta_data": meta_data},
                                 yield_kw_args)
    yield yield_class(**kw_args)


###############################################################################
#                    CONVERTING ALIGNMENTS TO MAF FORMAT                      #
###############################################################################

def sequence_to_maf_format(seq):
  res = ""

  if EMPTY_ALIGNMENT_STATUS_KEY in seq.meta_data:
    res += "e "
  else:
    res += "s "

  s = str(seq.start) if seq.is_positive_strand() else str(seq.remaining)
  res += (seq.name + " " + s + " " + str(seq.end - seq.start) + " " +
          seq.strand + " " + str(seq.remaining + seq.end) + " ")
  if EMPTY_ALIGNMENT_STATUS_KEY in seq.meta_data:
    res += seq.meta_data[EMPTY_ALIGNMENT_STATUS_KEY]
  else:
    res += seq.sequenceData
  if QUALITY_META_KEY in seq.meta_data:
    res += ("\n" + "q " + seq.name + " " + seq.meta_data[QUALITY_META_KEY])
  left_partial_present = (LEFT_STATUS_KEY in seq.meta_data or
                          LEFT_COUNT_KEY in seq.meta_data)
  right_partial_present = (RIGHT_STATUS_KEY in seq.meta_data or
                           RIGHT_COUNT_KEY in seq.meta_data)
  if left_partial_present or right_partial_present:
    # if one is present, all should be...
    left_partial_absent = (LEFT_STATUS_KEY not in seq.meta_data or
                           LEFT_COUNT_KEY not in seq.meta_data)
    right_partial_absent = (RIGHT_STATUS_KEY not in seq.meta_data or
                            RIGHT_COUNT_KEY not in seq.meta_data)
    if left_partial_absent or right_partial_absent:
      raise MAFError("converting alignment to MAF format string failed; " +
                     "found partially formed I-line information")
    res += ("\n" + "i " + seq.name + " ")
    res += (str(seq.meta_data[LEFT_STATUS_KEY]) + " " +
            str(seq.meta_data[LEFT_COUNT_KEY]) + " " +
            str(seq.meta_data[RIGHT_STATUS_KEY]) + " " +
            str(seq.meta_data[RIGHT_COUNT_KEY]))
  return res


def alignment_to_maf(alignment):
  res = ""

  # alignment line with meta data key-value pairs
  res += "a"
  for k in alignment.meta:
    if k == SEQ_ORDER_KEY:
      continue
    res += (" " + k + "=" + str(alignment.meta[k]))
  res += "\n"

  for k in alignment.meta[SEQ_ORDER_KEY]:
    res += (sequence_to_maf_format(alignment[k]) + "\n")
  return res


def write_maf(alig_iterable, output_stream):
  for alig in alig_iterable:
    output_stream.write(alignment_to_maf(alig))
    output_stream.write("\n")


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################

class TestMAF(unittest.TestCase):

  """Unit tests for MAF IO."""

  def setUp(self):
    """Set up some MAF files to use in unit tests."""
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
    self.b1_hg19 = Sequence("hg19.chr22", b1_hg19_seq, 1711, 1768,
                            "+", 51302798)
    self.b1_panTro = Sequence("panTro2.chrUn", b1_panTro_s, 1110, 1169, "+",
                              58616431 - 1169, {QUALITY_META_KEY: b1_panTro_q,
                                                LEFT_STATUS_KEY: "C",
                                                LEFT_COUNT_KEY: 0,
                                                RIGHT_STATUS_KEY: "C",
                                                RIGHT_COUNT_KEY: 0})
    self.b1_tarSyr = Sequence("tarSyr1.scaffold_5923", b1_tarSyr_s,
                              8928 - 2859 - 50, 8928 - 2859, "-",
                              2859, {QUALITY_META_KEY: b1_tarSyr_q,
                                     LEFT_STATUS_KEY: "N",
                                     LEFT_COUNT_KEY: 0,
                                     RIGHT_STATUS_KEY: "C",
                                     RIGHT_COUNT_KEY: 0})
    self.b1_mm4 = UnknownSequence("mm4.chr6", 53310102, 53310102 + 58, "+",
                                  151104725 - (53310102 + 58),
                                  {EMPTY_ALIGNMENT_STATUS_KEY: "I"})

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
    self.b2_hg19 = Sequence("hg19.chr22", b2_hg19_seq, 1772, 1825,
                            "+", 51302741)
    self.b2_panTro = Sequence("panTro2.chrUn", b2_panTro_s, 1169, 1169 + 53,
                              "+", 58616431 - (1169 + 53),
                              {QUALITY_META_KEY: b2_panTro_q,
                               LEFT_STATUS_KEY: "C",
                               LEFT_COUNT_KEY: 0,
                               RIGHT_STATUS_KEY: "C",
                               RIGHT_COUNT_KEY: 0})
    self.b2_tarSyr = Sequence("tarSyr1.scaffold_5923", b2_tarSyr_s,
                              8928 - 2909 - 124, 8928 - 2909, "-",
                              2909, {QUALITY_META_KEY: b2_tarSyr_q,
                                     LEFT_STATUS_KEY: "C",
                                     LEFT_COUNT_KEY: 0,
                                     RIGHT_STATUS_KEY: "N",
                                     RIGHT_COUNT_KEY: 0})
    self.b1 = b1
    self.b2 = b2

  def test_maf_iterator(self):
    """Test iterating over a MAF file."""
    blocks = [x for x in maf_iterator(StringIO.StringIO(self.maf1))]
    self.assertEqual(len(blocks), 2)
    self.assertEqual(blocks[0]["hg19.chr22"], self.b1_hg19)
    self.assertEqual(blocks[0]["panTro2.chrUn"], self.b1_panTro)
    self.assertEqual(blocks[0]["tarSyr1.scaffold_5923"], self.b1_tarSyr)
    self.assertEqual(blocks[0]["mm4.chr6"], self.b1_mm4)
    self.assertEqual(blocks[1]["hg19.chr22"], self.b2_hg19)
    self.assertEqual(blocks[1]["panTro2.chrUn"], self.b2_panTro)
    self.assertEqual(blocks[1]["tarSyr1.scaffold_5923"], self.b2_tarSyr)

    self.assertEqual(alignment_to_maf(blocks[0]).split(), self.b1.split())
    self.assertEqual(alignment_to_maf(blocks[1]).split(), self.b2.split())
    out = StringIO.StringIO()
    write_maf(maf_iterator(StringIO.StringIO(self.maf1)), out)
    self.assertEqual(self.maf1.split(), out.getvalue().split())

  def test_maf_comment_lines(self):
    """Simple test of MAF file that contains comment header line."""
    head = "##maf version=1 scoring=autoMZ.v1\n"
    blocks = [x for x in maf_iterator(StringIO.StringIO(head + self.maf1))]
    self.assertEqual(blocks[0]["hg19.chr22"], self.b1_hg19)


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
