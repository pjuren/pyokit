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


###############################################################################
#                                 CONSTANTS                                   #
###############################################################################

GAP_CHAR = '-'
UNKOWN_SEQ_NAME = "UNKNOWN_SEQUENCE"

# dafults values for reepatmasker formatting variables
DEFAULT_MAX_NAME_WIDTH = None
DEFAULT_COL_WIDTH = 50

# keys for meta-data dictionary; general
ANNOTATION_KEY = "annotation"
S1_NAME_KEY = "s1_name"
S2_NAME_KEY = "s2_name"
S1_START_KEY = "s1_start_coord"
S2_START_KEY = "s2_start_coord"
S1_END_KEY = "s1_end_coord"
S2_END_KEY = "s2_end_coord"
S1_START_NEG_STRAND_KEY = "s1_start_neg_coord"
S2_START_NEG_STRAND_KEY = "s2_start_neg_coord"
S1_END_NEG_STRAND_KEY = "s1_end_neg_coord"
S2_END_NEG_STRAND_KEY = "s2_end_neg_coord"
S1_REVERSE_COMP_KEY = "s1_rev_comp"
S2_REVERSE_COMP_KEY = "s2_rev_comp"
ALIG_SCORE_KEY = "alig_score"
PCENT_SUBS_KEY = "pcnt_subs"
PCENT_S1_INDELS_KEY = "pcnt_s1_indels"
PCENT_S2_INDELS_KEY = "pcnt_s2_indels"
# keys for meta-data dictionary; repeat-masker specific
UNKNOWN_RM_HEADER_FIELD_KEY = "unknown_rm_header_field"
RM_ID_KEY = "rm_id"
# keys for meta-data dictionary; special
ROUNDTRIP_KEY = "roundtrip"
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
    self._s1_ungapped_len = None

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
                         - ungapped_length(self.s1[i:i + column_width])
                         if self.s1_is_reverse_comp()
                         else s1_line_start_num
                         + ungapped_length(self.s1[i:i + column_width]))
      s2_line_start_num = (s2_line_end_num - 1 if self.s2_is_reverse_comp()
                           else s2_line_end_num + 1)
      s2_line_end_num = (s2_line_start_num
                         - ungapped_length(self.s2[i:i + column_width])
                         if self.s2_is_reverse_comp()
                         else s2_line_start_num
                         + ungapped_length(self.s2[i:i + column_width]))

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
