#!/usr/bin/python

"""
  Date of Creation: 11th Dec 2014

  Description:   Classes and functions for manipulating records of
                 retrotransposon occurrences

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
import unittest

# pyokit imports
from pyokit.datastruct.genomicInterval import GenomicInterval


###############################################################################
#                           MODULE-WIDE CONSTANTS                             #
###############################################################################

COMPLEMENT_CHARS = set(["C", "c"])


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class RetrotransposonError(Exception):
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


###############################################################################
#                 FUNCTIONS FOR PARSING STRING REPRESENTATIONS                #
###############################################################################

def from_repeat_masker_string(s):
  """
  Parse a RepeatMasker string. This format is white-space separated and has
  15 columns. The first 11 are invariably as follows:

  +----------------------+---------+-----------+----------------------------+
  |     Description      |  Type   |  Example  | Notes                      |
  +======================+=========+===========+============================+
  | alignment score      | int (?) |  463      |                            |
  +----------------------+---------+-----------+----------------------------+
  | percent divergence   | float   | 1.3       |                            |
  +----------------------+---------+-----------+----------------------------+
  | percent deletions    | float   | 0.6       |                            |
  +----------------------+---------+-----------+----------------------------+
  | percent insertions   | float   | 1.7       |                            |
  +----------------------+---------+-----------+----------------------------+
  | query seq. name      | string  | chr1      |                            |
  +----------------------+---------+-----------+----------------------------+
  | query seq. start     | int     | 10001     | always pos. strand coords  |
  +----------------------+---------+-----------+----------------------------+
  | query seq. end       | int     | 10468     | always pos. strand         |
  +----------------------+---------+-----------+----------------------------+
  | query seq. remaining | int     | (105)     | always in parenthesis,     |
  |                      |         |           | always negative strand     |
  |                      |         |           | coords. Bases left in      |
  |                      |         |           | query sequence before end  |
  |                      |         |           | (alt. can be interpreted   |
  |                      |         |           | as negative strand coord   |
  |                      |         |           | of genomic end coord).     |
  +----------------------+---------+-----------+----------------------------+
  | strand               | char    | C         | either + or C              |
  |                      |         |           | (complement); this refers  |
  |                      |         |           | to whether the match is    |
  |                      |         |           | to the consensus or the    |
  |                      |         |           | rev. comp. of it.          |
  +----------------------+---------+-----------+----------------------------+
  | repeat name          | string  | AluJo     |                            |
  +----------------------+---------+-----------+----------------------------+
  | repeat family        | string  | SINE/Alu  |                            |
  +----------------------+---------+-----------+----------------------------+

  The next three depend on whether the genomic sequence was a match to the
  consensus sequence, or to the reverse complement of the consensus sequence.
  In the former case, they are as follows:

  +------------------------+-------+-----------+-----------------------------+
  |     Description        | Type  |  Example  | Notes                       |
  +========================+=======+===========+=============================+
  | consensus match start  | int   |   45      | awlays pos strand           |
  +------------------------+-------+-----------+-----------------------------+
  | consensus match end    | int   |   60      | always pos strand;          |
  |                        |       |           | end >= start                |
  +------------------------+-------+-----------+-----------------------------+
  | consensus match remain | (int) |  (40)     | for '+' match, always       |
  |                        |       |           | negative strand coords. Num |
  |                        |       |           | of bases after match before |
  |                        |       |           | end of consensus. (Alt. can |
  |                        |       |           | be considered the negative  |
  |                        |       |           | strand coords of the match  |
  |                        |       |           | end in consensus seq.)      |
  +------------------------+-------+-----------+-----------------------------+

  In the latter case, they have this interpretation:

  +------------------------+-------+-----------+-----------------------------+
  |     Description        | Type  |  Example  | Notes                       |
  +========================+=======+===========+=============================+
  | consensus match before | (int) | (10)      | always neg. strand (always  |
  |                        |       |           | parens); number of bases    |
  |                        |       |           | before start of match on    |
  |                        |       |           | negative strand (alt. can   |
  |                        |       |           | be interpreted as the       |
  |                        |       |           | negative strand coordinates |
  |                        |       |           | of the match start)         |
  +------------------------+-------+-----------+-----------------------------+
  | consensus match start  | int   | 345       | start of match in pos.      |
  |                        |       |           | strand coords.              |
  +------------------------+-------+-----------+-----------------------------+
  | consensus match end    | int   | 143       | end of match in pos. strand |
  |                        |       |           | coords. Because match is on |
  |                        |       |           | neg. strand, 'start' is     |
  |                        |       |           | always >= end.              |
  +------------------------+-------+-----------+-----------------------------+

  The final column should be an ID, but I have noticed this is sometimes
  missing. We don't use all of the above fields, and we convert everything into
  positive strand coords where start < end always.

  :param s: string representation of a retrotransposon occurrence in
            RepeatMasker format.
  :return: Retrotransposon object corresponding to the string s.
  """
  parts = s.split()

  # get rid of those pesky parentheses
  for i in [11, 12, 13]:
    parts[i] = parts[i].strip()
    if parts[i][0] == "(":
      parts[i] = parts[i][1:]
    if parts[i][-1] == ")":
      parts[i] = parts[i][0:-1]

  if len(parts) != 15 and len(parts) != 14:
    raise RetrotransposonError("incorrectly formated RepeatMasker entry: " + s)

  # the first 11, invariant columns...
  # score = float(parts[0])  -- we don't use this
  # p_div = float(parts[1])  -- we don't use this
  # p_del = float(parts[2])  -- we don't use this
  # p_ins = float(parts[3])  -- we don't use this
  q_seq = parts[4].strip()
  q_seq_start = int(parts[5])
  q_seq_end = int(parts[6])
  # left_qseq = parts[7].strip() -- we don't use this
  strand = "-" if parts[8].strip() in COMPLEMENT_CHARS else "+"
  repeat_name = parts[9].strip()
  repeat_family = parts[10].strip()

  # the three columns the define match location in the consensus
  con_pos_start = int(parts[11]) if strand is "+" else int(parts[13])
  con_pos_end = int(parts[12])
  con_pos_left = int(parts[13]) if strand is "+" else int(parts[11])

  # might get ID, mgiht not -- its sometimes missing.
  # we won't be using this...
  # repeat_id = parts[14].strip() if len(parts) != 15 else None

  # we can compute the length of the consensus using the 'leftover' index
  con_seq_len = con_pos_end + con_pos_left

  rt = Retrotransposon(repeat_name, repeat_family)
  return RetrotransposonOccurrence(q_seq, q_seq_start, q_seq_end, strand,
                                   con_pos_start, con_pos_end, con_seq_len, rt)


###############################################################################
#                         THE Retrotransposon CLASS                           #
###############################################################################

class Retrotransposon(object):
  """
  Represents a retrotranspson; very simple, just name and family name. If
  avaialble, the consensus sequence for this RT.
  """

  def __init__(self, name, family_name, consensus=None):
    self.name = name
    self.family_name = family_name
    self.consensus = consensus


class RetrotransposonOccurrence(GenomicInterval):
  """
  Represents the occurrence of a retrotransposon in a genome.

  :param chrom:                  name of the sequence the occurrence is on
  :param genomic_start:          location in the sequence the occurrence starts
  :param genomic_end:            location in the sequence the occurrence ends
  :param consensus_match_strand: the strand of the consensus this occurrence
                                 is a match to ('+', or '-', where - means a
                                 match to the rev-comp of the consensus)
  :param consensus_start:
  :param consensus_end:
  :param concensus_len:
  :param retrotransposon:        the retrotransposon this is an occurrence of
  :param pairwise_alignment:     the full pairwise alignment between this
                                 occurrence and the consensus sequence
  """

  def __init__(self, chrom, genomic_start, genomic_end, consensus_match_strand,
               consensus_start, consensus_end, consensus_len, retrotransposon,
               pairwise_alignment=None):
    GenomicInterval.__init__(self, chrom, genomic_start, genomic_end,
                             retrotransposon.name, 0, "+")
    self.consensus_start = consensus_start
    self.consensus_end = consensus_end
    self.consensus_len = consensus_len
    self.consensus_match_strand = consensus_match_strand
    self.retrotransposon = retrotransposon
    # TODO would be good to provide some error-checking here to make sure the
    # alignment matches up properly with this occurrence...
    self.pairwise_alignment = pairwise_alignment

  def __str__(self):
    return GenomicInterval.__str__(self) + " --> " +\
        str(self.consensus_start) +\
        "\t" + str(self.consensus_end) + "\t" + str(self.consensus_len) +\
        "\t" + str(self.family_name) + "\t" +\
        str(self.consensus_match_strand)

  def liftover(self, intersecting_region):
    """
    Lift a region that overlaps the genomic occurrence of the retrotransposon
    to consensus sequence co-ordinates. If we ahve the full alignment, we will
    make use of that, otherwise we will fall back on the co-ordinates of the
    match start/end within the consensus. Note that since the alignment is
    gapped, the resultant co-ordinates can be slightly off if just using the
    start/end coordinates of the match

    :param intersecting_region: a region that intersects this occurrence.
    """

    # a little sanity check here to make sure intersecting_region really does..
    if not self.intersects(intersecting_region):
      raise RetrotransposonError("trying to lift " + str(intersecting_region) +
                                 " from genomic to transposon coordinates " +
                                 "in " + str(self) + ", but it doesn't " +
                                 "intersect!")

    if self.pariwise_alignment is not None:
      s, e = self.pairwise_alignment.liftover(intersecting_region.start,
                                              intersecting_region.end)
    else:
      s = max(intersecting_region.start - self.start, 0)
      e = min(s + (intersecting_region.end - intersecting_region.start),
              self.consensus_len)
      if self.consensus_match_strand is "-":
        e = self.consensus_len - s
        s = max(e - (intersecting_region.end - intersecting_region.start),
                0)

    return GenomicInterval(self.name, s, e,
                           intersecting_region.name, intersecting_region.score,
                           self.strand)


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestRetrotransposon(unittest.TestCase):
  """
  Test cases for this module.
  """

  def test_from_rm_string(self):
    """
    test the conversion of repeat-masker strings to retrotransposon occurrence
    objects; note that a round-trip test is not possible because we don't store
    enough of the RM fields to reconstruct it.
    """

    rm_1 = ("463   1.3  0.6  1.7  chr1        10001   10468 (249240153) +  "
            + "(TAACCC)n      Simple_repeat            1  463    (0)      1")
    rm_2 = ("3612  11.4 21.5  1.3  chr14        10469   11447 (249239174) C  "
            + "TAR1           Satellite/telo       (399) 1712    483      2")
    rm_3 = ("484  25.1 13.2  0.0  chrX        11505   11675 (249238946) C  "
            + "L1MC5a         LINE/L1             (2382) 5648   5452      3")

    rm_1_p = from_repeat_masker_string(rm_1)
    rm_2_p = from_repeat_masker_string(rm_2)
    rm_3_p = from_repeat_masker_string(rm_3)

    self.assertEqual(rm_1_p.chrom, "chr1")
    self.assertEqual(rm_2_p.chrom, "chr14")
    self.assertEqual(rm_3_p.chrom, "chrX")

    self.assertEqual(rm_1_p.start, 10001)
    self.assertEqual(rm_2_p.start, 10469)
    self.assertEqual(rm_3_p.start, 11505)

    self.assertEqual(rm_1_p.end, 10468)
    self.assertEqual(rm_2_p.end, 11447)
    self.assertEqual(rm_3_p.end, 11675)

    self.assertEqual(rm_1_p.consensus_match_strand, '+')
    self.assertEqual(rm_2_p.consensus_match_strand, '-')
    self.assertEqual(rm_3_p.consensus_match_strand, '-')

    self.assertEqual(rm_1_p.consensus_start, 1)
    self.assertEqual(rm_2_p.consensus_start, 483)
    self.assertEqual(rm_3_p.consensus_start, 5452)

    self.assertEqual(rm_1_p.consensus_end, 463)
    self.assertEqual(rm_2_p.consensus_end, 1712)
    self.assertEqual(rm_3_p.consensus_end, 5648)

    self.assertEqual(rm_1_p.consensus_len, 463)
    self.assertEqual(rm_2_p.consensus_len, 2111)
    self.assertEqual(rm_3_p.consensus_len, 8030)

    self.assertEqual(rm_1_p.retrotransposon.name, "(TAACCC)n")
    self.assertEqual(rm_2_p.retrotransposon.name, "TAR1")
    self.assertEqual(rm_3_p.retrotransposon.name, "L1MC5a")

    self.assertEqual(rm_1_p.retrotransposon.family_name, "Simple_repeat")
    self.assertEqual(rm_2_p.retrotransposon.family_name, "Satellite/telo")
    self.assertEqual(rm_3_p.retrotransposon.family_name, "LINE/L1")

    self.assertEqual(rm_1_p.retrotransposon.consensus, None)
    self.assertEqual(rm_2_p.retrotransposon.consensus, None)
    self.assertEqual(rm_3_p.retrotransposon.consensus, None)

    self.assertEqual(rm_1_p.pairwise_alignment, None)
    self.assertEqual(rm_2_p.pairwise_alignment, None)
    self.assertEqual(rm_3_p.pairwise_alignment, None)


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == '__main__':
    unittest.main()
