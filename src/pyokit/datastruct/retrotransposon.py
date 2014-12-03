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

  alignment score      --  int (?) -- e.g.: 463
  percent divergence   --  float   -- e.g.: 1.3
  percent deletions    --  float   -- e.g.: 0.6
  percent insertions   --  float   -- e.g.: 1.7
  query seq. name      --  string  -- e.g.: chr1
  query seq. start     --  int     -- e.g.: 10001  -- always pos. strand coords
  query seq. end       --  int     -- e.g.: 10468  -- always pos. strand
  query seq. remaining --  (int)   -- e.g.: (105)  -- always in parenthesis,
                                                      always negative strand
                                                      coords. Bases left in
                                                      query sequence before end
                                                      (alt. can be interpreted
                                                      as negative strand coord
                                                      of genomic end coord).
  strand               --  char    -- e.g.: +      -- either + or C
                                                      (complement); this refers
                                                      to whether the match is
                                                      to the consensus or the
                                                      rev. comp. of it.
  repeat name          --  string  -- e.g. AluJo
  repeat family        --  string  -- e.g. SINE/Alu

  The next three depend on whether the genomic sequence was a match to the
  consensus sequence, or to the reverse complement of the consensus sequence.
  In the former case, they are as follows:

  consensus match start  -- int   -- e.g.: 45   -- awlays pos strand;
  consensus match end    -- int   -- e.g.: 60   -- always pos strand;
                                                   end >= start
  consensus match remain -- (int) -- e.g.: (40) -- for '+' match, always
                                                   negative strand coords. Num
                                                   of bases after match before
                                                   end of consensus. (Alt. can
                                                   be considered the negative
                                                   strand coords of the match
                                                   end in consensus seq.)

  In the latter case, they have this interpretation:

  consensus match before -- (int) -- e.g.: (10) -- always neg. strand (always
                                                   parens); number of bases
                                                   before start of match on
                                                   negative strand (alt. can
                                                   be interpreted as the
                                                   negative strand coordinates
                                                   of the match start)
  consensus match start  -- int   -- e.g.: 345  -- start of match in pos.
                                                   strand coords.
  consensus match end    -- int   -- e.g.: 143  -- end of match in pos. strand
                                                   coords. Because match is on
                                                   neg. strand, 'start' is
                                                   always >= end.

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

  return Retrotransposon(q_seq, q_seq_start, q_seq_end, strand, con_pos_start,
                         con_pos_end, con_seq_len, repeat_name, repeat_family)


###############################################################################
#                         THE Retrotransposon CLASS                           #
###############################################################################

class Retrotransposon(GenomicInterval):
  """
  Represents the occurrence of a retrotransposon in a genome.

  :param genomic_interval:
  :param consensus_start:
  :param consensus_end:
  :param name:
  :param family_name:
  """

  def __init__(self, chrom, genomic_start, genomic_end, consensus_match_strand,
               consensus_start, consensus_end, consensus_len, name,
               family_name):
    GenomicInterval.__init__(self, chrom, genomic_start, genomic_end,
                             name, 0, "+")
    self.consensus_start = consensus_start
    self.consensus_end = consensus_end
    self.consensus_len = consensus_len
    self.family_name = family_name
    self.consensus_match_strand = consensus_match_strand

  def __str__(self):
    return GenomicInterval.__str__(self) + " --> " +\
        str(self.consensus_start) +\
        "\t" + str(self.consensus_end) + "\t" + str(self.consensus_len) +\
        "\t" + str(self.family_name) + "\t" +\
        str(self.consensus_match_strand)

  def liftover(self, intersecting_region):
    """
    Lift a region that overlaps the genomic occurrence of the retrotransposon
    to consensus sequence co-ordinates. Note that since the alignment is
    gapped, the resultant co-ordinates can be slightly off. This could be
    fixed if we had the full details of the alignment, but I haven't
    implemented that yet... maybe later.
    """

    # a little sanity check here to make sure intersecting_region really does..
    if not self.intersects(intersecting_region):
      raise RetrotransposonError("trying to lift " + str(intersecting_region) +
                                 " from genomic to transposon coordinates " +
                                 "in " + str(self) + ", but it doesn't " +
                                 "intersect!")

    start = max(intersecting_region.start - self.start, 0)
    end = min(start + (intersecting_region.end - intersecting_region.start),
              self.consensus_len)
    if self.consensus_match_strand is "-":
      end = self.consensus_len - start
      start = max(end - (intersecting_region.end - intersecting_region.start),
                  0)

    return GenomicInterval(self.name, start, end,
                           intersecting_region.name, intersecting_region.score,
                           self.strand)
