"""
  Date of Creation: 11th Dec 2014

  Description:   Classes and for representing pairwise and multiple sequence
                 alignments.

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

###############################################################################
#                                 CONSTANTS                                   #
###############################################################################

GAP_CHAR = '-'
UNKOWN_SEQ_NAME = "UNKNOWN_SEQUENCE"

# dafults values for reepatmasker formatting variables
DEFAULT_MAX_NAME_WIDTH = None
DEFAULT_COL_WIDTH = 50


###############################################################################
#                   KEYS FOR META-DATA DICTIONARY; GENERAL                    #
###############################################################################

#: name of the first sequence
S1_NAME_KEY = "s1_name"

#: name of the second sequence
S2_NAME_KEY = "s2_name"

#: start co-ordinate of the first sequence
S1_START_KEY = "s1_start_coord"

#: start co-ordinate of the second sequence
S2_START_KEY = "s2_start_coord"

#: end co-ordinate of the first sequence
S1_END_KEY = "s1_end_coord"

#: end co-ordinate of the second sequence
S2_END_KEY = "s2_end_coord"

#: start co-ordinate of the first sequenc in negative strand co-ordinates
S1_START_NEG_STRAND_KEY = "s1_start_neg_coord"

#: start co-ordinate of the second sequenc in negative strand co-ordinates
S2_START_NEG_STRAND_KEY = "s2_start_neg_coord"

#: end co-ordinate of the first sequenc in negative strand co-ordinates
S1_END_NEG_STRAND_KEY = "s1_end_neg_coord"

#: end co-ordinate of the second sequenc in negative strand co-ordinates
S2_END_NEG_STRAND_KEY = "s2_end_neg_coord"

#: boolean -- is the first sequence reverse complemented?
S1_REVERSE_COMP_KEY = "s1_rev_comp"

#: boolean -- is the second sequence reverse complemented?
S2_REVERSE_COMP_KEY = "s2_rev_comp"

#: score for the alignment; no requirement on how this is computed.
ALIG_SCORE_KEY = "alig_score"

#: annotation for each column of the alignment; no requirements on format.
ANNOTATION_KEY = "annotation"

#: percentage of substitutions (i.e. columns that have no gaps, not matches)
PCENT_SUBS_KEY = "pcnt_subs"

#: the percentage of the first sequence which is gaps. This is not computed
#: from the sequence itself, so may not be accurate.
PCENT_S1_INDELS_KEY = "pcnt_s1_indels"

#: the percentage of the second sequence which is gaps. This is not computed
#: from the sequence itself, so may not be accurate.
PCENT_S2_INDELS_KEY = "pcnt_s2_indels"

#:
ROUNDTRIP_KEY = "roundtrip"

###############################################################################
#          KEYS FOR META-DATA DICTIONARY; SEPCIFIC TO REPEAT-MASKER           #
###############################################################################

#: If this is a repeat-masker alignment, this stores the value of the
#: second-last token from the header of the alignment; its meaning is unknown.
UNKNOWN_RM_HEADER_FIELD_KEY = "unknown_rm_header_field"

#: If this is a repeat-masker alignment, this stores the unique ID assigned
#: to it by repeat-masker
RM_ID_KEY = "rm_id"

###############################################################################
#                               FULL KEY LIST                                 #
###############################################################################

#: The full set of meta-data keys that have special meaning
KNOWN_KEYS = set([ANNOTATION_KEY, S1_NAME_KEY, S1_END_KEY, S2_END_KEY,
                  S1_START_NEG_STRAND_KEY, S2_START_NEG_STRAND_KEY,
                  S1_END_NEG_STRAND_KEY, S2_END_NEG_STRAND_KEY, S2_NAME_KEY,
                  S1_START_KEY, S2_START_KEY, S1_REVERSE_COMP_KEY,
                  S2_REVERSE_COMP_KEY, RM_ID_KEY, UNKNOWN_RM_HEADER_FIELD_KEY,
                  PCENT_SUBS_KEY, PCENT_S1_INDELS_KEY, PCENT_S2_INDELS_KEY,
                  ALIG_SCORE_KEY])


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class AlignmentError(Exception):
  """
  Class representing errors that occur when manipulating pairwise or multiple
  alignment objects
  """
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


###############################################################################
#               PAIRWISE ALIGNMENT CLASS AND HELPER FUNCTIONS                 #
###############################################################################


def ungapped_length(s):
  """
  compute the length of a sequence without counting gap characters.

  :param s: the sequence to compute the lenght of
  :return: integer, the number of non-gap characters in the sequence
  """
  res = 0
  for c in s:
    if c != GAP_CHAR:
      res += 1
  return res


class PairwiseAlignment(object):
  """
  An alignment of two sequences (DNA, RNA, protein...).

  :param s1: the first sequence, with gaps
  :param s2: the second sequence, with gaps
  :param meta_data: a dictionary with key-value pairs representing meta-data
                    about this alignment.
  """

  def __init__(self, s1, s2, meta_data=None):
    self.s1 = s1
    self.s2 = s2
    if len(s1) != len(s2):
      raise AlignmentError("Invalid pairwise alignment, sequences have "
                           "different lengths: '" + str(s1) + "' -and- '" +
                           str(s2) + "'")
    self.meta = meta_data

    # handling names; if we don't have sequence names, then just set the meta
    # data values to empty strings
    if S1_NAME_KEY not in self.meta:
      self.meta[S1_NAME_KEY] = ""
    if S2_NAME_KEY not in self.meta:
      self.meta[S2_NAME_KEY] = ""

    # handling coordinates; if we don't know coordinates for the sequences,
    # set these to start at position 1
    if S1_START_KEY not in self.meta:
      self.meta[S1_START_KEY] = 1
    if S2_START_KEY not in self.meta:
      self.meta[S2_START_KEY] = 1

    # we will compute ungapped lengths of the sequences on-demand
    self._s1_ungapped_len = None
    self._s2_ungapped_len = None

  @property
  def s1_ungapped_len(self):
    """
    The length of the first sequence in this alignment, not counting gaps
    """
    if self._s1_ungapped_len is None:
      self._s1_ungapped_len = ungapped_length(self.s1)
    return self._s1_ungapped_len

  @property
  def s2_ungapped_len(self):
    """
    The length of the second sequence in this alignment, not counting gaps
    """
    if self._s2_ungapped_len is None:
      self._s2_ungapped_len = ungapped_length(self.s2)
    return self._s2_ungapped_len

  def s1_is_reverse_comp(self):
    """
    :return: True if the first sequence in this alignment is on the negative
             strand, False otherwise
    """
    return S1_REVERSE_COMP_KEY in self.meta and self.meta[S1_REVERSE_COMP_KEY]

  def s2_is_reverse_comp(self):
    """
    :return: True if the second sequence in this alignment is on the negative
             strand, False otherwise
    """
    return S2_REVERSE_COMP_KEY in self.meta and self.meta[S2_REVERSE_COMP_KEY]

  def size(self):
    """
    return the size of this alignment (defined as the number of columns in the
    alignemnt).
    """
    return len(self.s1)

  def repeat_masker_header(self):
    """
    generate the header string of a repeatmasker formated representation of
    this pairwise alignment.
    """

    res = ""
    res += str(self.meta[ALIG_SCORE_KEY]) + " "
    res += "{:.2f}".format(self.meta[PCENT_SUBS_KEY]) + " "
    res += "{:.2f}".format(self.meta[PCENT_S1_INDELS_KEY]) + " "
    res += "{:.2f}".format(self.meta[PCENT_S2_INDELS_KEY]) + " "
    res += (self.meta[S1_NAME_KEY] if S1_NAME_KEY in self.meta
            else UNKOWN_SEQ_NAME) + " "
    res += str(self.meta[S1_START_KEY]) + " "
    res += str(self.meta[S1_END_KEY]) + " "
    res += "(" + str(self.meta[S1_END_NEG_STRAND_KEY]) + ") "
    res += ("C " if self.s2_is_reverse_comp() else "")
    res += (self.meta[S2_NAME_KEY] if S2_NAME_KEY in self.meta
            else UNKOWN_SEQ_NAME) + " "
    res += ("(" + str(self.meta[S2_START_NEG_STRAND_KEY]) + ")"
            if self.s2_is_reverse_comp() else str(self.meta[S2_START_KEY]))
    res += " "
    res += (str(self.meta[S2_START_KEY]) if self.s2_is_reverse_comp()
            else str(self.meta[S2_END_KEY])) + " "
    res += (str(self.meta[S2_END_KEY]) if self.s2_is_reverse_comp()
            else "(" + str(self.meta[S2_END_NEG_STRAND_KEY]) + ")") + " "
    res += self.meta[UNKNOWN_RM_HEADER_FIELD_KEY] + " "
    res += self.meta[RM_ID_KEY]
    return res

  def __str__(self):
    """
    return a string representation of this pairwise alignment
    """
    return self.to_repeat_masker_string()

  def sequence_to_alignment_coords(self, seq_num, start, end):
    """
    convert an interval in one of the sequences into an interval in the
    alignment. Alignment intervals are inclusive of start, but not end. They
    are zero-based. Hence the full alignment has coords [0, N), where N is the
    length of the alignment (number of columns). Sequence coords follow the
    same conventions: zero-based, inclusive of start but not end.

    :param seq_num: which sequence are the start and end coords for? 1 or 2
    :param start:   start of the interval in sequence co-ordinates
    :param end:     end of the interval in sequence co-ordinates
    """
    assert(seq_num == 1 or seq_num == 2)
    assert(end > start)
    try:
      s_start = (self.meta[S1_START_KEY] if seq_num == 1
                 else self.meta[S2_START_KEY])
      s_end = (self.meta[S1_END_KEY] if seq_num == 1
               else self.meta[S2_END_KEY])
      seq = self.s1 if seq_num == 1 else self.s2
    except KeyError:
      s_start = 0
      s_end = self.s1_ungapped_length - 1
    assert(s_start <= start < s_end)
    assert(s_start < end <= s_end)

    num_gaps = 0
    num_non_gaps = 0
    res = []
    current_start = None
    current_end = None
    print seq
    for i in range(0, end):
      if seq[i] == GAP_CHAR:
        num_gaps += 1
      else:
        num_non_gaps += 1

      if (num_non_gaps > end - s_start):
        # we're done, gone past the end of the ROI
        break
      if (num_non_gaps > start - s_start):
        # within ROI
        print "in ROI at position " + str(i) + " (" + seq[i] + ")"
        if seq[i] != GAP_CHAR:
          if current_start == None and current_end == None:
            current_start = i
            current_end = i + 1
          else:
            if seq[i - 1] == GAP_CHAR:
              # is the start of a new non-gapped region...
              print "\tstarts new non-gapped region"
              print "\t\tadd (" + str(current_start) + "," + str(current_end) + ")"
              res.append((current_start, current_end))
              current_start = i
              current_end = i + 1
            if seq[i - 1] != GAP_CHAR:
              # is continuation of non-gapped region
              print "\tcontinue non-gapped region"
              current_end += 1
    res.append((current_start, current_end))
    return res

  def to_repeat_masker_string(self, column_width=DEFAULT_COL_WIDTH,
                              m_name_width=DEFAULT_MAX_NAME_WIDTH):
    """
    generate a repeatmasker formated representation of this pairwise alignment.

    :param column_width: number of characters to output per line of alignment
    :param m_name_width: truncate names on alignment lines to this length
                         (set to None for no truncation)
    """
    # figure out the complement column
    s1_comp = "C" if self.s1_is_reverse_comp() else " "
    s2_comp = "C" if self.s2_is_reverse_comp() else " "

    # figure out the maximum name length, so we can size that column properly;
    # pre-compute the space-padded names too
    s1_len = len(self.meta[S1_NAME_KEY])
    s2_len = len(self.meta[S2_NAME_KEY])
    f_len = max(s1_len, s2_len)
    if m_name_width != None:
      f_len = min(f_len, m_name_width)
    s1_n = self.meta[S1_NAME_KEY][:f_len] + (' ' * (f_len - s1_len))
    s2_n = self.meta[S2_NAME_KEY][:f_len] + (' ' * (f_len - s2_len))

    # figure out the max width for the coordinates column; we use size of the
    # alignment here rather than ungapped coordinates because its an upper
    # bound and easier to compute (i.e. for sure already know).
    s1_line_end_num = (self.meta[S1_START_KEY] + 1 if self.s1_is_reverse_comp()
                       else self.meta[S1_START_KEY] - 1)
    s2_line_end_num = (self.meta[S2_START_KEY] + 1 if self.s2_is_reverse_comp()
                       else self.meta[S2_START_KEY] - 1)
    max_num_len = max(len(str(self.meta[S1_START_KEY] + self.size())),
                      len(str(self.meta[S2_START_KEY] + self.size())))

    res = ""  # our result
    i = 0     # how much of the full, gapped alignment, has been output so far?
    res += self.repeat_masker_header() + "\n\n"
    while i < len(self.s1):
      # keep track of how much of each sequence we've output
      s1_line_start_num = (s1_line_end_num - 1 if self.s1_is_reverse_comp()
                           else s1_line_end_num + 1)
      s1_line_end_num = (s1_line_start_num
                         - ungapped_length(self.s1[i:i + column_width]) + 1
                         if self.s1_is_reverse_comp()
                         else s1_line_start_num
                         + ungapped_length(self.s1[i:i + column_width]) - 1)
      s2_line_start_num = (s2_line_end_num - 1 if self.s2_is_reverse_comp()
                           else s2_line_end_num + 1)
      s2_line_end_num = (s2_line_start_num
                         - ungapped_length(self.s2[i:i + column_width]) + 1
                         if self.s2_is_reverse_comp()
                         else s2_line_start_num
                         + ungapped_length(self.s2[i:i + column_width]) - 1)

      # output sequence one
      res += (s1_comp + " " + s1_n + " ")
      s1_line_start_num_str = str(s1_line_start_num)
      s1_num_padding = max_num_len - len(s1_line_start_num_str)
      res += (' ' * s1_num_padding) + s1_line_start_num_str + " "
      res += self.s1[i:i + column_width] + " "
      res += str(s1_line_end_num) + "\n"

      # output the annotation string, if we have one; needs to be padded by the
      # number of char in the name col (f_len), the number in the coordinate
      # col (max_num_len), the one char in the complement columns, and the
      # three spaces that are used as column seperators for those.
      if ANNOTATION_KEY in self.meta :
        res += (((f_len + max_num_len) * ' ') + "    " +
                self.meta[ANNOTATION_KEY][i:i + column_width] + "\n")

      # output sequence two
      res += (s2_comp + " " + s2_n + " ")
      s2_line_start_num_str = str(s2_line_start_num)
      s2_num_padding = max_num_len - len(s2_line_start_num_str)
      res += (' ' * s2_num_padding) + s2_line_start_num_str + " "
      res += self.s2[i:i + column_width] + " "
      res += str(s2_line_end_num) + "\n"

      res += "\n"
      i += column_width

    # otuput any meta data key-value pairs that aren't known to us.
    if self.meta != None:
      for k in self.meta:
        if k not in KNOWN_KEYS:
          if k is ROUNDTRIP_KEY:
            res += (self.meta[k] + "\n")
          else:
            res += (k + " = " + str(self.meta[k]) + "\n")

    # remove any trailing whitespace
    res = res.strip()
    return res


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################
class TestAlignments(unittest.TestCase):

  def test_sequence_to_alig_coord(self):
    """
    test converting co-ordinates for an interval within a component sequence
    of a pairwise alignment into co-ordinates within the alignment itself
    """
    meta = {}
    meta[S1_START_KEY] = 100
    meta[S2_START_KEY] = 1000
    meta[S1_END_KEY] = 129
    meta[S2_END_KEY] = 934
    meta[S1_REVERSE_COMP_KEY] = False
    meta[S2_REVERSE_COMP_KEY] = True

    p = PairwiseAlignment("-TCGCGTAGC---CGC-TAGCTGATGCGAT-CTGA",
                          "ATCGCGTAGCTAGCGCG-AGCTG---CGATGCT--", meta)

    self.failUnlessEqual(p.sequence_to_alignment_coords(1, 103, 111),
                         [(4, 10), (13, 15)])


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
