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
import math

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
  |                        |       |           | end >= start; it is         |
  |                        |       |           | inclusive of the final      |
  |                        |       |           | coordinate (unlike e.g.     |
  |                        |       |           | UCSC coords)                |
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

  if len(parts) != 15 and len(parts) != 14:
    raise RetrotransposonError("incorrectly formated RepeatMasker entry: " + s)

  # get rid of those pesky parentheses
  for i in [11, 12, 13]:
    parts[i] = parts[i].strip()
    if parts[i][0] == "(":
      parts[i] = parts[i][1:]
    if parts[i][-1] == ")":
      parts[i] = parts[i][0:-1]

  # the first 11, invariant columns...
  # score = float(parts[0])  -- we don't use this
  # p_div = float(parts[1])  -- we don't use this
  # p_del = float(parts[2])  -- we don't use this
  # p_ins = float(parts[3])  -- we don't use this
  q_seq = parts[4].strip()
  q_seq_start = int(parts[5])
  q_seq_end = int(parts[6]) + 1
  # left_qseq = parts[7].strip() -- we don't use this
  strand = "-" if parts[8].strip() in COMPLEMENT_CHARS else "+"
  repeat_name = parts[9].strip()
  repeat_family = parts[10].strip()

  # the three columns the define match location in the consensus
  con_pos_start = int(parts[11]) if strand is "+" else int(parts[13])
  con_pos_end = int(parts[12]) + 1
  con_pos_left = int(parts[13]) if strand is "+" else int(parts[11])

  # might get ID, mgiht not -- its sometimes missing.
  repeat_id = int(parts[14].strip()) if len(parts) == 15 else None

  # we can compute the length of the consensus using the 'leftover' index
  con_seq_len = con_pos_end + con_pos_left

  rt = Retrotransposon(repeat_name, repeat_family)
  return RetrotransposonOccurrence(q_seq, q_seq_start, q_seq_end, strand,
                                   con_pos_start, con_pos_end, con_seq_len,
                                   rt, uniq_id=repeat_id)


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
  Represents the occurrence of a retrotransposon in a genome. Note that
  occurrences are always consdiered to be on the positive strand of the
  genome, and might match either the concensus sequence of the repeat, or
  its reverse complement.

  :param chrom:                  name of the sequence the occurrence is on
  :param genomic_start:          location in the sequence the occurrence starts
  :param genomic_end:            location in the sequence the occurrence ends
  :param consensus_match_strand: the strand of the consensus this occurrence
                                 is a match to ('+', or '-', where - means a
                                 match to the rev-comp of the consensus)
  :param consensus_start:        start coordinate of the match, in consensus
                                 sequence coordinates.
  :param consensus_end:          end coordinate of the match, in consensus
                                 sequence coordinates.
  :param concensus_len:          the length of the consensus sequence (not
                                 the length of the match in the consensus
                                 sequence)
  :param retrotransposon:        the retrotransposon this is an occurrence of
  :param pairwise_alignment:     the full pairwise alignment between this
                                 occurrence and the consensus sequence
  :param uniq_id:                a unique identifier for this alignment, if
                                 one exists
  """

  def __init__(self, chrom, genomic_start, genomic_end, consensus_match_strand,
               consensus_start, consensus_end, consensus_len, retrotransposon,
               pairwise_alignment=None, uniq_id=None):
    """
    Constructor for RetrotransposonOccurrence objects; see class docstring
    for parameter defs.
    """
    GenomicInterval.__init__(self, chrom, genomic_start, genomic_end,
                             retrotransposon.name, 0, "+")
    self.consensus_start = consensus_start
    self.consensus_end = consensus_end
    # TODO check that consensus length is at least as long as end index
    self.consensus_len = consensus_len
    self.consensus_match_strand = consensus_match_strand
    self.retrotransposon = retrotransposon
    self.uniq_id = uniq_id
    # TODO would be good to provide some error-checking here to make sure the
    # alignment matches up properly with this occurrence...
    self.pairwise_alignment = pairwise_alignment

  def __str__(self):
    """
    :return: a string representation of this occurrence in the format
             <genomic_interval_details> --> match_start match_end \
                                            consensus_len repeat_name \
                                            match_strand
    """
    return GenomicInterval.__str__(self) + " --> " +\
        str(self.consensus_start) + "\t" + str(self.consensus_end) + "\t" +\
        str(self.consensus_len) + "\t" + self.repeat_name() + "\t" +\
        str(self.consensus_match_strand)

  def repeat_name(self):
    """
    :return: the name of the repeat that this is an occurrence in, using the
             format family#sub-family
    """
    return (self.retrotransposon.family_name + "#" + self.retrotransposon.name)

  def liftover(self, intersecting_region):
    """
    Lift a region that overlaps the genomic occurrence of the retrotransposon
    to consensus sequence co-ordinates. This method will behave differently
    depending on whether this retrotransposon occurrance contains a full
    alignment or not. If it does, the alignment is used to do the liftover and
    an exact result is provided. If it does not, the coordinates are used to
    do the liftover, padding either the genomic region or consensus sequence
    (whichever is shorter) with equally spaced gaps to make the size of both
    match.

    :param intersecting_region: a region that intersects this occurrence.
    :return: list of GenomicInterval objects. This is a list because a genomic
             deletion of part of the retrotransposon can fragment the
             intersecting region and result in more than one returned interval.
    """

    # a little sanity check here to make sure intersecting_region really does..
    if not self.intersects(intersecting_region):
      raise RetrotransposonError("trying to lift " + str(intersecting_region) +
                                 " from genomic to transposon coordinates " +
                                 "in " + str(self) + ", but it doesn't " +
                                 "intersect!")

    if self.pairwise_alignment is not None:
      return self.pairwise_alignment.liftover(self.chrom, self.repeat_name(),
                                              intersecting_region.start,
                                              intersecting_region.end,
                                              trim=True)
    return self.liftover_coordinates(intersecting_region)

  def liftover_coordinates(self, intersecting_region):
    """
    Lift a region that overlaps the genomic occurrence of the retrotransposon
    to consensus sequence co-ordinates. Use only the coordinates of the
    occurrence and the consensus match (i.e. ignore a full alignment, even
    if one is present). If the genomic size of the match is the same as the
    size of the match in the consensus sequence then it is assumed no indels
    are present in either sequence, otherwise indels are assumed to exist
    uniformly distributed along the sequences to make their lengths match.

    :param intersect_region: a region that intersects this occurrence.
    :return: list of GenomicInterval objects. This is a list because a genomic
             deletion of part of the retrotransposon can fragment the
             intersecting region and result in more than one returned interval.
    """
    # a little sanity check here to make sure intersecting_region really does..
    if not self.intersects(intersecting_region):
      raise RetrotransposonError("trying to lift " + str(intersecting_region) +
                                 " from genomic to transposon coordinates " +
                                 "in " + str(self) + ", but it doesn't " +
                                 "intersect!")

    consensus_match_length = self.consensus_end - self.consensus_start
    size_dif = consensus_match_length - len(self)

    if size_dif == 0:
      # put this in a list, even though we know we're getting only one
      # interval, to make the return type consistent
      return [self.__liftover_coordinates_size_match(intersecting_region)]
    else:
      return self.__liftover_coordinates_genomic_indels(intersecting_region)

  def __liftover_coordinates_size_match(self, intersecting_region):
    """
    Lift a region that overlaps the genomic occurrence of the retrotransposon
    to consensus sequence coordinates using just the coordinates (not the full)
    alignment, when they have the same length. This is an internal helper
    method. The above length constraint must be true, otherwise an assertion
    is failed.

    :param intersect_region: a region that intersects this occurrence.
    :return: GenomicInterval representing the region after lifting to consensus
    :note: no checks are made for whether the interval really intersects!!
    """
    # should always pass, but check anyway...
    consensus_match_length = self.consensus_end - self.consensus_start
    assert(consensus_match_length - len(self) == 0)

    if self.consensus_match_strand is '+':
      s = max(intersecting_region.start - self.start, 0) +\
          self.consensus_start
      e = min(max(intersecting_region.end - self.start, 0) +
              self.consensus_start, self.consensus_len)
      g = GenomicInterval(self.repeat_name(), s, e, intersecting_region.name,
                          intersecting_region.score, self.strand)
      return g
    elif self.consensus_match_strand is '-':
      e = (self.consensus_end -
           max(intersecting_region.start - self.start, 0))
      s = (self.consensus_end -
           min(max(intersecting_region.end - self.start, 0), len(self)))
      g = GenomicInterval(self.repeat_name(), s, e, intersecting_region.name,
                          intersecting_region.score, self.strand)
      return g
    else:
      raise RetrotransposonError("couldn't determine strand of " +
                                 "retrotransposon occurrance " + str(self))

  def __liftover_coordinates_genomic_indels(self, intersect_region):
    """
    Lift a region that overlaps the genomic occurrence of the retrotransposon
    to consensus sequence coordinates using just the coordinates (not the full
    alignment), when they have differing length. This is an internal helper
    method. The above length constraint must be true, otherwise an assertion
    is failed.

    :param intersect_region: a region that intersects this occurrence.
    :return: list of GenomicInterval objects. This is a list because a genomic
             deletion of part of the retrotransposon can fragment the
             intersecting region and result in more than one returned interval.
             List might be empty, if intersection overlaps only gaps.
    :note: no checks are made for whether the interval really intersects!!
    """
    # should never happen, but check anyway...
    consensus_match_length = self.consensus_end - self.consensus_start
    size_dif = consensus_match_length - len(self)
    assert(size_dif != 0)

    if size_dif < 0:
      # genomic region longer; assume genomic region contains insertions
      return self.__liftover_coordinates_genomic_insertions(intersect_region)
    else:
      # genomic region shorter; assume genomic region contains deletions
      return self.__liftover_coordinates_genomic_deletions(intersect_region)

  def __liftover_coordinates_genomic_insertions(self, intersecting_region):
    """
    Lift a region that overlaps the genomic occurrence of this repeat to the
    consensus sequence coordinates using just coordinates (not the full
    alignment, even if it is avaialble), when the length of the genomic match
    is greater than the concensus match. We assume there are insertions in
    the genomic sequence (i.e. gaps in the consensus match). We uniformly
    distribute the gaps through the consensus match. e.g. if the genomic
    region is 100nt long and the consensus match is 90nt long, we will consider
    an insertion to occurr every 10nt in the genomic sequence. The start and
    end coordinates of the region after lifting will be reduced by the number
    of genomic insertions that would have come before them.

    :param intersecting_region: a region that intersects this occurrence.
    :return: list of GenomicInterval objects. This is a list because, although
             the liftover can produce only one genomic interval in consensus,
             it might also produce none (if the interval overlaps only
             insertions in the genomic sequence), hence an empty list is
             possible.
    """
    # should never happen, but check anyway...
    consensus_match_length = self.consensus_end - self.consensus_start
    size_dif = consensus_match_length - len(self)
    assert(size_dif < 0)

    gap_interval = len(self) / (-1 * size_dif)
    s_dist_to_gen_start = max(intersecting_region.start - self.start, 0)
    e_dist_to_gen_start = max(intersecting_region.end - self.start, 0)

    if self.consensus_match_strand is '+':
      s = s_dist_to_gen_start + self.consensus_start
      e = e_dist_to_gen_start + self.consensus_start
      s = s - (s_dist_to_gen_start / gap_interval)
      e = min(e - (e_dist_to_gen_start / gap_interval), self.consensus_len)
    else:
      e = self.consensus_end - s_dist_to_gen_start
      s = self.consensus_end - e_dist_to_gen_start
      s = max(s + (e_dist_to_gen_start / gap_interval),
              self.consensus_start)
      e = min(e + (s_dist_to_gen_start / gap_interval), self.consensus_end)

    res = [] if s == e else [GenomicInterval(self.repeat_name(), s, e,
                                             intersecting_region.name,
                                             intersecting_region.score,
                                             self.strand)]
    return res

  def __liftover_coordinates_genomic_deletions(self, intersecting_region):
    """
    A 'private' helper member function to perform liftover in coordinate space
    when the length of the genomic match is smaller than the concensus match.
    We assume the genomic region contains deletions. In this case, we uniformly
    distribute the deletions (gaps) through the genomic region. This is the
    trickiest liftover case because it fragments the region after liftover.
    This method should only be called when the aboe condition is true (longer
    consensus region than genomic region) otherwise an assertion will be
    failed.

    :param intersecting_region: a region that intersects this occurrence.
    :return: list of GenomicInterval objects. This is a list because a genomic
             deletion of part of the retrotransposon can fragment the
             intersecting region and result in more than one returned interval.
    """
    # should never happen, but check anyway...
    consensus_match_length = self.consensus_end - self.consensus_start
    size_dif = consensus_match_length - len(self)
    assert(size_dif > 0)

    name = self.repeat_name()
    gap_interval = int(math.ceil(len(self) / float(size_dif)))
    s_dist_to_gen_start = max(intersecting_region.start - self.start, 0)
    e_dist_to_gen_start = min(max(intersecting_region.end -
                                  self.start, 0), len(self))
    gaps_before = s_dist_to_gen_start / gap_interval
    gaps_in = ((e_dist_to_gen_start - 1) / gap_interval) - gaps_before
    gen_s_dist = s_dist_to_gen_start
    left_to_lift = (min(self.end, intersecting_region.end) -
                    max(self.start, intersecting_region.start))
    res = []
    if self.consensus_match_strand is '+':
      s = s_dist_to_gen_start + self.consensus_start
      s = s + gaps_before
      for i in range(0, gaps_in):
        e = s + min((gap_interval - (gen_s_dist % gap_interval)), left_to_lift)
        res.append(GenomicInterval(name, s, e, intersecting_region.name,
                                   intersecting_region.score, self.strand))
        gen_s_dist += (e - s)
        left_to_lift -= (e - s)
        s = e + 1
      e = min(s + min(gap_interval, left_to_lift), self.consensus_end)
      if e - s != 0:
        res.append(GenomicInterval(name, s, e, intersecting_region.name,
                                   intersecting_region.score, self.strand))
    else:
      e = self.consensus_end - s_dist_to_gen_start
      e = e - gaps_before
      for i in range(0, gaps_in):
        s = e - min((gap_interval - (gen_s_dist % gap_interval)), left_to_lift)
        res.append(GenomicInterval(name, s, e, intersecting_region.name,
                                   intersecting_region.score, self.strand))
        gen_s_dist += (e - s)
        left_to_lift -= (e - s)
        e = s - 1
      s = max(e - min(gap_interval, left_to_lift), self.consensus_start)
      if e - s != 0:
        res.append(GenomicInterval(name, s, e, intersecting_region.name,
                                   intersecting_region.score, self.strand))
    return res


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestRetrotransposon(unittest.TestCase):
  """
  Test cases for this module.
  """

  def setUp(self):
    self.rt1 = Retrotransposon("(TAACCC)n", "Simple_repeat")
    self.rt2 = Retrotransposon("TAR1", "Satellite/telo")
    self.rt3 = Retrotransposon("L1MC5a", "LINE/L1")
    self.rto0 = RetrotransposonOccurrence("chr1", 10001, 10468, '+', 1, 468,
                                          468, self.rt1, uniq_id=0)
    self.rto1 = RetrotransposonOccurrence("chr1", 10001, 10468, '+', 1, 463,
                                          463, self.rt1, uniq_id=1)
    self.rto2 = RetrotransposonOccurrence("chr14", 10469, 11447, '-', 483,
                                          1712, 1229, self.rt2, uniq_id=2)
    self.rto3 = RetrotransposonOccurrence("chrX", 11505, 11675, '-',
                                          5452, 5648, 196, self.rt3, uniq_id=3)
    self.rto4 = RetrotransposonOccurrence("chr14", 10469, 11447, '-', 469,
                                          1447, 3000, self.rt2, uniq_id=4)
    self.rto5 = RetrotransposonOccurrence("chrX", 11505, 11675, '-',
                                          505, 701, 1000, self.rt3, uniq_id=5)
    self.rto6 = RetrotransposonOccurrence("chr5", 10, 50, '-',
                                          10, 40, 196, self.rt3, uniq_id=6)
    self.rto7 = RetrotransposonOccurrence("chr5", 10, 18, '+',
                                          110, 121, 196, self.rt3, uniq_id=7)

    self.rm_0 = ("463   1.3  0.6  1.7  chr1      10001   10467 (249240153) +" +
                 "  (TAACCC)n      Simple_repeat          1  467    (0)     1")
    self.rm_1 = ("463   1.3  0.6  1.7  chr1      10001   10467 (249240153) +" +
                 "  (TAACCC)n      Simple_repeat          1  462    (0)     1")
    self.rm_2 = ("3612  11.4 21.5  1.3  chr14    10469   11446 (249239174) C" +
                 "  TAR1           Satellite/telo     (399) 1711    483     2")
    self.rm_3 = ("484  25.1 13.2  0.0  chrX      11505   11674 (249238946) C" +
                 "  L1MC5a         LINE/L1           (2382) 5647   5452     3")
    self.rm_4 = ("3612  11.4 21.5  1.3  chr14    10469   11447 (249239174) C" +
                 "  TAR1           Satellite/telo     (399) 1447    467     2")
    self.rm_5 = ("484  25.1 13.2  0.0  chrX      11505   11675 (249238946) C" +
                 "  L1MC5a         LINE/L1           (2382) 675   503     3")
    self.rm_6 = ("146  25.6 14.9  0.0  chr5      10 49 (249238946) C " +
                 "L1MC5a         LINE/L1             (2382) 10    39      6")

  def test_from_rm_string(self):
    """
    test the conversion of repeat-masker strings to retrotransposon occurrence
    objects; note that a round-trip test is not possible because we don't store
    enough of the RM fields to reconstruct it.
    """

    rm_1_p = from_repeat_masker_string(self.rm_1)
    rm_2_p = from_repeat_masker_string(self.rm_2)
    rm_3_p = from_repeat_masker_string(self.rm_3)

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

  def test_liftover_pos_strand_no_alignment_no_indels(self):
    """
    test lifting regions from genomic coordinates into repeat coordinates when
    the match is to the consensus sequence (not the reverse complement) and
    we don't have a full allignment to use

    the test case that can be used for this are:
    rto0   "chr1",  10001   10468  '+'   ->    1  468      +

    (others are all matches to reverse complement of consensus); the first
    gives use a mapping where the size of the genomic region is the same as
    the matching region in the consensus, while the second gives us a mapping
    where the genomic region is larger than the consensus and so presumably
    contains some insertions, and finally the last one gives us a mapping
    where the genomic region is ...
    """

    # intersecting region is fully contained within repeat.
    in1 = GenomicInterval("chr1", 10005, 10015, "AA", 50, '+')
    lf1 = self.rto0.liftover(in1)
    expc1 = [GenomicInterval("Simple_repeat#(TAACCC)n", 5, 15, "AA", 50, '+')]
    self.assertEqual(lf1, expc1)

    # intersecting region overlaps the start of the match.
    in2 = GenomicInterval("chr1", 9995, 10015, "BB", 10, '+')
    lf2 = self.rto0.liftover(in2)
    expc2 = [GenomicInterval("Simple_repeat#(TAACCC)n", 1, 15, "BB", 10, '+')]
    self.assertEqual(lf2, expc2)

    # we have no full alignment, the match is the same size as the consensus
    # and the intersecting region is at the start of the match (boundary case)
    in3 = GenomicInterval("chr1", 10001, 10015, "CC", 20, '+')
    lf3 = self.rto0.liftover(in3)
    expc3 = [GenomicInterval("Simple_repeat#(TAACCC)n", 1, 15, "CC", 20, '+')]
    self.assertEqual(lf3, expc3)

    # intersecting region overlaps the end of the match
    in4 = GenomicInterval("chr1", 10460, 10470, "DD", 14, '+')
    lf4 = self.rto0.liftover(in4)
    expc4 = GenomicInterval("Simple_repeat#(TAACCC)n", 460, 468, "DD", 14, '+')
    self.assertEqual(lf4, [expc4])

    # intersecting region is at the end of the match (boundary case)
    in5 = GenomicInterval("chr1", 10467, 10470, "EE", 50, '+')
    lf5 = self.rto0.liftover(in5)
    expc5 = GenomicInterval("Simple_repeat#(TAACCC)n", 467, 468, "EE", 50, '+')
    self.assertEqual(lf5, [expc5])

    # we have no full alignment, the match is the same size as the consensus
    # and the intersecting region ends before (at) the match; this should fail
    # with an exception
    in6 = GenomicInterval("chr1", 9995, 10001, "FF", 50, '+')
    self.assertRaises(RetrotransposonError, self.rto0.liftover, in6)

    # we have no full alignment, the match is the same size as the consensus
    # and the intersecting region start after (at) the match; this should fail
    # with an exception
    in7 = GenomicInterval("chr1", 10468, 10470, "GG", 50, '+')
    self.assertRaises(RetrotransposonError, self.rto0.liftover, in7)

  def test_liftover_neg_strand_no_alignment_no_indels(self):
    """
    test lifting regions from genomic coordinates to repeat coordinates when
    the repeat occurrence is a match to the reverse complement of the
    consensus.

    rto4 chr14 10469 11447  -->  469  1447  -

    """
    # intersecting region is fully contained within repeat.
    in1 = GenomicInterval("chr14", 10470, 10480, "AA", 50, '+')
    lf1 = self.rto4.liftover(in1)
    expc1 = [GenomicInterval("Satellite/telo#TAR1", 1436, 1446, "AA",
                             50, '+')]
    self.assertEqual(lf1, expc1)

    # intersecting region overlaps the start of the match.
    in2 = GenomicInterval("chr14", 10465, 10475, "BB", 10, '+')
    lf2 = self.rto4.liftover(in2)
    expc2 = [GenomicInterval("Satellite/telo#TAR1", 1441, 1447, "BB",
                             10, '+')]
    self.assertEqual(lf2, expc2)

    # intersecting region is at the start of the match (boundary case)
    in3 = GenomicInterval("chr14", 10459, 10470, "CC", 20, '+')
    lf3 = self.rto4.liftover(in3)
    expc3 = [GenomicInterval("Satellite/telo#TAR1", 1446, 1447, "CC", 20, '+')]
    self.assertEqual(lf3, expc3)

    # intersecting region overlaps the end of the match
    in4 = GenomicInterval("chr14", 11440, 11450, "DD", 14, '+')
    lf4 = self.rto4.liftover(in4)
    expc4 = GenomicInterval("Satellite/telo#TAR1", 469, 476, "DD", 14, '+')
    self.assertEqual(lf4, [expc4])

    # intersecting region is at the end of the match (boundary case)
    in5 = GenomicInterval("chr14", 11446, 11450, "EE", 50, '+')
    lf5 = self.rto4.liftover(in5)
    expc5 = GenomicInterval("Satellite/telo#TAR1", 469, 470, "EE", 50, '+')
    self.assertEqual(lf5, [expc5])

    # intersecting region ends before (at) the match; this should fail
    # with an exception
    in6 = GenomicInterval("chr14", 9995, 10469, "FF", 50, '+')
    self.assertRaises(RetrotransposonError, self.rto4.liftover, in6)

    # intersecting region start after (at) the match; this should fail
    # with an exception
    in7 = GenomicInterval("chr14", 11447, 11470, "GG", 50, '+')
    self.assertRaises(RetrotransposonError, self.rto4.liftover, in7)

  def test_liftover_no_alignment_genomic_insertions(self):
    """
    test lifting regions when there is an unknown net insertion in the genomic
    region (i.e. genomic region is larger than match to consensus)

    rto1 -->   chr1  10001  10468  -->  1, 463, +
    [5 bp larger; space = 467 / 5 = 93]
    """

    # intersecting region overlaps start and no insertions in genomic seq.
    in1 = GenomicInterval("chr1", 9990, 10010, "AA", 50, '+')
    lf1 = self.rto1.liftover(in1)
    expc1 = [GenomicInterval("Simple_repeat#(TAACCC)n", 1, 10, "AA", 50, '+')]
    self.assertEqual(lf1, expc1)

    # intersecting region overlaps start and some insertions in genomic seq.
    in2 = GenomicInterval("chr1", 9990, 10200, "BB", 50, '+')
    lf2 = self.rto1.liftover(in2)
    expc2 = [GenomicInterval("Simple_repeat#(TAACCC)n", 1, 198, "BB", 50, '+')]
    self.assertEqual(lf2, expc2)

    # intersecting region overlaps first insertion in genomic sequence
    in3 = GenomicInterval("chr1", 10090, 10100, "CC", 50, '+')
    lf3 = self.rto1.liftover(in3)
    expc3 = [GenomicInterval("Simple_repeat#(TAACCC)n", 90, 99, "CC", 50, '+')]
    self.assertEqual(lf3, expc3)

    # intersecting region overlaps whole match
    in4 = GenomicInterval("chr1", 9990, 10500, "DD", 50, '+')
    lf4 = self.rto1.liftover(in4)
    expc4 = [GenomicInterval("Simple_repeat#(TAACCC)n", 1, 463, "DD", 50, '+')]
    self.assertEqual(lf4, expc4)

    # intersecting region overlaps end and some insertions
    in5 = GenomicInterval("chr1", 10460, 10470, "EE", 50, '+')
    lf5 = self.rto1.liftover(in5)
    expc5 = [GenomicInterval("Simple_repeat#(TAACCC)n", 456, 463,
                             "EE", 50, '+')]
    self.assertEqual(lf5, expc5)

    # intersecting region overlaps end and no insertions
    in6 = GenomicInterval("chr1", 10466, 10470, "FF", 50, '+')
    lf6 = self.rto1.liftover(in6)
    expc6 = [GenomicInterval("Simple_repeat#(TAACCC)n", 461, 463,
                             "FF", 50, '+')]
    self.assertEqual(lf6, expc6)

    # start sits on a gap
    in7 = GenomicInterval("chr1", 10093, 10096, "GG", 50, '+')
    lf7 = self.rto1.liftover(in7)
    expc7 = [GenomicInterval("Simple_repeat#(TAACCC)n", 93, 95,
                             "GG", 50, '+')]
    self.assertEqual(lf7, expc7)

    # start is first after gap
    in8 = GenomicInterval("chr1", 10094, 10096, "GG", 50, '+')
    lf8 = self.rto1.liftover(in8)
    expc8 = [GenomicInterval("Simple_repeat#(TAACCC)n", 93, 95,
                             "GG", 50, '+')]
    self.assertEqual(lf8, expc8)

    # end sits on a gap
    in9 = GenomicInterval("chr1", 10090, 10093, "GG", 50, '+')
    lf9 = self.rto1.liftover(in9)
    expc9 = [GenomicInterval("Simple_repeat#(TAACCC)n", 90, 93,
                             "GG", 50, '+')]
    self.assertEqual(lf9, expc9)

    # region covers only gaps -- should get an empty list
    in10 = GenomicInterval("chr1", 10093, 10094, "FF", 50, '+')
    lf10 = self.rto1.liftover(in10)
    expc10 = []
    self.assertEqual(lf10, expc10)

  def test_liftover_no_alignment_genomic_insertions_neg(self):
    """
    test lifting regions when there is an unknown net insertion in the genomic
    region (i.e. genomic region is larger than match to consensus) and the
    match is to the reverse complement of the consensus sequence

    rto6  chr5 10 50 --> 10 40 -   (LINE/L1#L1MC5a)
    [40 to 30; genomic region is larger by 10bp; insert space is 4]
    """

    # intersecting region overlaps start and no insertions in genomic seq.
    in1 = GenomicInterval("chr5", 8, 12, "AA", 50, '+')
    lf1 = self.rto6.liftover(in1)
    expc1 = [GenomicInterval("LINE/L1#L1MC5a", 38, 40, "AA", 50, '+')]
    self.assertEqual(lf1, expc1)

    # intersecting region overlaps start and some insertions in genomic seq.
    in2 = GenomicInterval("chr5", 8, 20, "BB", 50, '+')
    lf2 = self.rto6.liftover(in2)
    expc2 = [GenomicInterval("LINE/L1#L1MC5a", 32, 40, "BB", 50, '+')]
    self.assertEqual(lf2, expc2)

    # intersecting region overlaps whole match
    in3 = GenomicInterval("chr5", 8, 55, "CC", 50, '+')
    lf3 = self.rto6.liftover(in3)
    expc3 = [GenomicInterval("LINE/L1#L1MC5a", 10, 40, "CC", 50, '+')]
    self.assertEqual(lf3, expc3)

    # intersecting region overlaps end and some insertions
    in4 = GenomicInterval("chr5", 40, 55, "DD", 50, '+')
    lf4 = self.rto6.liftover(in4)
    expc4 = [GenomicInterval("LINE/L1#L1MC5a", 10, 17, "DD", 50, '+')]
    self.assertEqual(lf4, expc4)

    # start sits on a gap
    in5 = GenomicInterval("chr5", 14, 17, "EE", 50, '+')
    lf5 = self.rto6.liftover(in5)
    expc5 = [GenomicInterval("LINE/L1#L1MC5a", 34, 37, "EE", 50, '+')]
    self.assertEqual(lf5, expc5)

    # start is first after gap
    in6 = GenomicInterval("chr5", 15, 17, "FF", 50, '+')
    lf6 = self.rto6.liftover(in6)
    expc6 = [GenomicInterval("LINE/L1#L1MC5a", 34, 36, "FF", 50, '+')]
    self.assertEqual(lf6, expc6)

    # start is first before gap
    in7 = GenomicInterval("chr5", 13, 17, "GG", 50, '+')
    lf7 = self.rto6.liftover(in7)
    expc7 = [GenomicInterval("LINE/L1#L1MC5a", 34, 37, "GG", 50, '+')]
    self.assertEqual(lf7, expc7)

    # end is first before gap
    in8 = GenomicInterval("chr5", 10, 13, "HH", 50, '+')
    lf8 = self.rto6.liftover(in8)
    expc8 = [GenomicInterval("LINE/L1#L1MC5a", 37, 40, "HH", 50, '+')]
    self.assertEqual(lf8, expc8)

    # end sits on gap
    in9 = GenomicInterval("chr5", 10, 14, "II", 50, '+')
    lf9 = self.rto6.liftover(in9)
    expc9 = [GenomicInterval("LINE/L1#L1MC5a", 37, 40, "II", 50, '+')]
    self.assertEqual(lf9, expc9)

    # end is first after gap
    in10 = GenomicInterval("chr5", 10, 15, "JJ", 50, '+')
    lf10 = self.rto6.liftover(in10)
    expc10 = [GenomicInterval("LINE/L1#L1MC5a", 36, 40, "JJ", 50, '+')]
    self.assertEqual(lf10, expc10)

    # region covers only gaps -- should get an empty list
    in11 = GenomicInterval("chr5", 13, 14, "KK", 50, '+')
    lf11 = self.rto6.liftover(in11)
    self.assertEqual(lf11, [])

  def test_liftover_no_alignment_genomic_deletions(self):
    """
    test lifting regions when there is an unknown net deletion in the genomic
    region (i.e. genomic region is shorter than match to consensus)

    rto7  chr5 10  18  -->  110  121  +
    [8 to 11; genomic region smaller by 3bp -- gap interval is 3]
    10  11  12  -   13  14  15  -   16  17  -
    110 111 112 113 114 115 116 117 118 119 120
    """

    # intersecting region is fully contained in match, covers some gaps;
    # starts partway through a link, finishes on link boundary
    in1 = GenomicInterval("chr5", 11, 14, "AA", 50, '+')
    lf1 = self.rto7.liftover(in1)
    expc1 = [GenomicInterval("LINE/L1#L1MC5a", 111, 113, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 114, 115, "AA", 50, '+')]
    self.assertEqual(lf1, expc1)

    # intersecting region fully covers the match
    in2 = GenomicInterval("chr5", 8, 20, "AA", 50, '+')
    lf2 = self.rto7.liftover(in2)
    expc2 = [GenomicInterval("LINE/L1#L1MC5a", 110, 113, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 114, 117, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 118, 120, "AA", 50, '+')]
    self.assertEqual(lf2, expc2)

    # intersecting region overlaps start and some gaps; finished partway
    # through a link
    in3 = GenomicInterval("chr5", 8, 14, "AA", 50, '+')
    lf3 = self.rto7.liftover(in3)
    expc3 = [GenomicInterval("LINE/L1#L1MC5a", 110, 113, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 114, 115, "AA", 50, '+')]
    self.assertEqual(lf3, expc3)

    # intersecting region overlaps end and some gaps
    in4 = GenomicInterval("chr5", 13, 20, "AA", 50, '+')
    lf4 = self.rto7.liftover(in4)
    expc4 = [GenomicInterval("LINE/L1#L1MC5a", 114, 117, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 118, 120, "AA", 50, '+')]
    self.assertEqual(lf4, expc4)

  def test_liftover_no_alignment_genomic_deletions_neg(self):
    """
    test lifting regions when there is an unknown net deletion in the genomic
    region (i.e. genomic region is shorter than match to consensus) and the
    match is to the reverse complement of the consensus sequence

    rto3  chrX   11505  11675  -->  5452 5648  -
    [170 to 196; genomic region is smaller by 26bp -- gap interval is 7]
    """

    # intersecting region is fully contained in match, covers some gaps
    in1 = GenomicInterval("chrX", 11508, 11520, "AA", 50, '+')
    lf1 = self.rto3.liftover(in1)
    expc1 = [GenomicInterval("LINE/L1#L1MC5a", 5641, 5645, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5633, 5640, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5631, 5632, "AA", 50, '+')]
    self.assertEqual(lf1, expc1)

    # intersecting region fully covers the match
    in2 = GenomicInterval("chrX", 11500, 11700, "AA", 50, '+')
    lf2 = self.rto3.liftover(in2)
    expc2 = [GenomicInterval("LINE/L1#L1MC5a", 5641, 5648, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5633, 5640, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5625, 5632, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5617, 5624, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5609, 5616, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5601, 5608, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5593, 5600, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5585, 5592, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5577, 5584, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5569, 5576, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5561, 5568, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5553, 5560, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5545, 5552, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5537, 5544, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5529, 5536, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5521, 5528, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5513, 5520, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5505, 5512, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5497, 5504, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5489, 5496, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5481, 5488, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5473, 5480, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5465, 5472, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5457, 5464, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5454, 5456, "AA", 50, '+')]
    self.assertEqual(lf2, expc2)

    # intersecting region overlaps start and some gaps; finished partway
    # through a link
    in3 = GenomicInterval("chrX", 11500, 11515, "AA", 50, '+')
    lf3 = self.rto3.liftover(in3)
    expc3 = [GenomicInterval("LINE/L1#L1MC5a", 5641, 5648, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5637, 5640, "AA", 50, '+')]
    self.assertEqual(lf3, expc3)

    # intersecting region overlaps end and some gaps
    in4 = GenomicInterval("chrX", 11665, 11700, "AA", 50, '+')
    lf4 = self.rto3.liftover(in4)
    expc4 = [GenomicInterval("LINE/L1#L1MC5a", 5465, 5466, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5457, 5464, "AA", 50, '+'),
             GenomicInterval("LINE/L1#L1MC5a", 5454, 5456, "AA", 50, '+')]
    self.assertEqual(lf4, expc4)


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == '__main__':
    unittest.main()
