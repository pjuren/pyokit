"""
  Date of Creation: 11th Dec 2014

  Description:   Classes and functions for IO of repeatmasker pairwise
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
import sys
import StringIO
import unittest
import re

# Pyokit imports
from pyokit.datastruct.sequence import Sequence
from pyokit.datastruct.sequence import UNKNOWN_SEQ_NAME
from pyokit.util.progressIndicator import ProgressIndicator
from pyokit.datastruct.multipleAlignment import PairwiseAlignment
from pyokit.datastruct.multipleAlignment import JustInTimePairwiseAlignment


###############################################################################
#                                 CONSTANTS                                   #
###############################################################################

REPEATMASKER_FIELDS_TO_TRIM = [0, 1, 4, 7, 8, 9, 11]
REPEATMASKER_VALIDATE_MUTATIONS = True
REPEATMASKER_VALID_ANN_CHARS = ['-', 'i', 'v', ' ', '?']

# dafults values for reepatmasker formatting variables
DEFAULT_MAX_NAME_WIDTH = None
DEFAULT_COL_WIDTH = 50


###############################################################################
#                        META-DATA ALIGNMENT CONSTANTS                        #
#  these are the keys to index into a pairwise alignment object's meta-data   #
#  dictionary and extract whatever is needed for formatting a repeat-masker   #
#                            representation of it                             #
###############################################################################

# stores the value of the second-last token from the header of the alignment;
# its meaning is unknown to me.
UNKNOWN_RM_HEADER_FIELD_KEY = "unknown_rm_header_field"

# stores the unique ID assigned to an alignment but repeat-masker
RM_ID_KEY = "rm_id"

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
#                               FULL KEY LIST                                 #
###############################################################################

#: The full set of meta-data keys that have special meaning
KNOWN_KEYS = set([ANNOTATION_KEY, PCENT_SUBS_KEY, PCENT_S1_INDELS_KEY,
                  PCENT_S2_INDELS_KEY, ALIG_SCORE_KEY, RM_ID_KEY,
                  UNKNOWN_RM_HEADER_FIELD_KEY])


###############################################################################
#           CONVERTING ALIGNMENTS TO REPEATMASKER STRING FORMATS              #
###############################################################################

def _get_repeat_masker_header(pairwise_alignment):
  """generate header string of repeatmasker formated repr of self."""
  res = ""
  res += str(pairwise_alignment.meta[ALIG_SCORE_KEY]) + " "
  res += "{:.2f}".format(pairwise_alignment.meta[PCENT_SUBS_KEY]) + " "
  res += "{:.2f}".format(pairwise_alignment.meta[PCENT_S1_INDELS_KEY]) + " "
  res += "{:.2f}".format(pairwise_alignment.meta[PCENT_S2_INDELS_KEY]) + " "
  res += (pairwise_alignment.s1.name
          if (pairwise_alignment.s1.name != "" and
              pairwise_alignment.s1.name is not None)
          else UNKNOWN_SEQ_NAME) + " "
  res += str(pairwise_alignment.s1.start) + " "
  res += str(pairwise_alignment.s1.end - 1) + " "
  res += "(" + str(pairwise_alignment.s1.remaining) + ") "
  res += ("C " if not pairwise_alignment.s2.is_positive_strand() else "")
  res += (pairwise_alignment.s2.name
          if (pairwise_alignment.s2.name != "" and
              pairwise_alignment.s2.name is not None)
          else UNKNOWN_SEQ_NAME) + " "
  res += ("(" + str(pairwise_alignment.s2.remaining) + ")"
          if not pairwise_alignment.s2.is_positive_strand()
          else str(pairwise_alignment.s2.start))
  res += " "
  # Note here that we need to convert between our internal representation
  # for coordinates and the repeat-masker one; internally, we always store
  # coordinates as exclusive of the final value with start < end;
  # repeatmasker gives the larger coordinate as the 'start' when the match
  # is to the reverse complement, so we have to swap start/end, and its
  # coordinates are inclusive of end, so we have to subtract 1 from end.
  res += str(pairwise_alignment.s2.end - 1) + " "
  res += (str(pairwise_alignment.s2.start)
          if not pairwise_alignment.s2.is_positive_strand()
          else "(" + str(pairwise_alignment.s2.remaining) + ")") + " "
  res += pairwise_alignment.meta[UNKNOWN_RM_HEADER_FIELD_KEY] + " "
  res += str(pairwise_alignment.meta[RM_ID_KEY])
  return res


def _to_repeatmasker_string(pairwise_alignment, column_width=DEFAULT_COL_WIDTH,
                            m_name_width=DEFAULT_MAX_NAME_WIDTH):
  """
  generate a repeatmasker formated representation of this pairwise alignment.

  :param column_width: number of characters to output per line of alignment
  :param m_name_width: truncate names on alignment lines to this length
                       (set to None for no truncation)
  """
  s1 = pairwise_alignment.s1
  s2 = pairwise_alignment.s2

  s1_neg = not s1.is_positive_strand()
  s2_neg = not s2.is_positive_strand()
  size = pairwise_alignment.size()

  # figure out the complement column
  s1_comp = "C" if s1_neg else " "
  s2_comp = "C" if s2_neg else " "

  # figure out the maximum name length, so we can size that column properly;
  # pre-compute the space-padded names too
  s1_len = len(s1.name)
  s2_len = len(s2.name)
  f_len = max(s1_len, s2_len)
  if m_name_width is not None:
    f_len = min(f_len, m_name_width)
  s1_n = s1.name[:f_len] + (' ' * (f_len - s1_len))
  s2_n = s2.name[:f_len] + (' ' * (f_len - s2_len))

  # figure out the max width for the coordinates column; we use size of the
  # alignment here rather than ungapped coordinates because its an upper
  # bound and easier to compute (i.e. for sure already know).
  s1_line_end_num = (s1.end if s1_neg else s1.start - 1)
  s2_line_end_num = (s2.end if s2_neg else s2.start - 1)
  max_num_len = max(len(str(s1.start + size)), len(str(s2.start + size)))

  res = ""  # our result
  i = 0     # how much of the full, gapped alignment, has been output so far?
  res += _get_repeat_masker_header(pairwise_alignment) + "\n\n"
  while i < len(pairwise_alignment.s1):
    # keep track of how much of each sequence we've output
    s1_sub = s1.gapped_relative_subsequence(i + 1, min(i + column_width + 1, len(s1) + 1))
    s2_sub = s2.gapped_relative_subsequence(i + 1, min(i + column_width + 1, len(s2) + 1))
    s1_ug_len = s1_sub.ungapped_len
    s2_ug_len = s2_sub.ungapped_len
    s1_line_start_num = (s1_line_end_num - 1 if s1_neg
                         else s1_line_end_num + 1)
    s1_line_end_num = (s1_line_start_num - s1_ug_len + 1 if s1_neg
                       else s1_line_start_num + s1_ug_len - 1)
    s2_line_start_num = (s2_line_end_num - 1 if s2_neg
                         else s2_line_end_num + 1)
    s2_line_end_num = (s2_line_start_num - s2_ug_len + 1 if s2_neg
                       else s2_line_start_num + s2_ug_len - 1)

    # output sequence one
    res += (s1_comp + " " + s1_n + " ")
    s1_line_start_num_str = str(s1_line_start_num)
    s1_num_padding = max_num_len - len(s1_line_start_num_str)
    res += (' ' * s1_num_padding) + s1_line_start_num_str + " "
    res += pairwise_alignment.s1[i:i + column_width] + " "
    res += str(s1_line_end_num) + "\n"

    # output the annotation string, if we have one; needs to be padded by the
    # number of char in the name col (f_len), the number in the coordinate
    # col (max_num_len), the one char in the complement columns, and the
    # three spaces that are used as column seperators for those.
    if ANNOTATION_KEY in pairwise_alignment.meta:
      res += (((f_len + max_num_len) * ' ') + "    " +
              pairwise_alignment.meta[ANNOTATION_KEY][i:i + column_width] + "\n")

    # output sequence two
    res += (s2_comp + " " + s2_n + " ")
    s2_line_start_num_str = str(s2_line_start_num)
    s2_num_padding = max_num_len - len(s2_line_start_num_str)
    res += (' ' * s2_num_padding) + s2_line_start_num_str + " "
    res += pairwise_alignment.s2[i:i + column_width] + " "
    res += str(s2_line_end_num) + "\n"

    res += "\n"
    i += column_width

  # otuput any meta data key-value pairs that aren't known to us.
  if pairwise_alignment.meta is not None:
    for k in pairwise_alignment.meta:
      if k not in KNOWN_KEYS:
        if k is ROUNDTRIP_KEY:
          res += (pairwise_alignment.meta[k] + "\n")
        else:
          res += (k + " = " + str(pairwise_alignment.meta[k]) + "\n")

  # remove any trailing whitespace
  res = res.strip()
  return res


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class AlignmentIteratorError(Exception):
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


###############################################################################
#         REPEATMASKER PAIRWISE ALIGNMENT ITERATOR AND HELPER FUNCS           #
###############################################################################

def _rm_is_alignment_line(parts, s1_name, s2_name):
  """
  return true if the tokenized line is a repeatmasker alignment line.

  :param parts:   the line, already split into tokens around whitespace
  :param s1_name: the name of the first sequence, as extracted from the header
                  of the element this line is in
  :param s2_name: the name of the second sequence, as extracted from the header
                  of the element this line is in
  """
  if len(parts) < 2:
    return False
  if _rm_name_match(parts[0], s1_name):
    return True
  if (_rm_name_match(parts[0], s2_name) or
     (parts[0] == "C" and _rm_name_match(parts[1], s2_name))):
    return True
  return False


def _rm_is_header_line(parts, n):
  """
  determine whether a pre-split string is a repeat-masker alignment header.

  headers have no special structure or symbol to mark them, so this is based
  only on the number of elements, and what data type they are.
  """
  if (n == 15 and parts[8] == "C"):
    return True
  if (n == 14 and parts[0].isdigit()):
    return True


def _rm_is_valid_annotation_line(line):
  """
  :return: True if the line contains only valid annotation characters (defined
           in REPEATMASKER_VALID_ANN_CHARS), otherwise False
  """
  for c in line:
    if c not in REPEATMASKER_VALID_ANN_CHARS:
      return False
  return True


def _rm_compute_leading_space_alig(space_pres_split, seq):
  """
  count the number of characters that precede the sequence in a repeatmasker
  alignment line. E.g. in the following line:
         '  chr1               11 CCCTGGAGATTCTTATT--AGTGATTTGGGCT 41'
  the answer would be 24.

  :param space_pres_split: the alignment line, split into tokens around spaces,
                           but with the spaces conserved as tokens.
  :param seq: the sequence token.
  """
  c = 0
  for i in range(0, len(space_pres_split)):
    if space_pres_split[i] == seq:
      break
    c += len(space_pres_split[i])
  return c


def _rm_compute_leading_space(space_s_pres_split):
  """
  count the number of spaces that precede a non-space token (not including
  empty string tokens) in a string.

  :param space_s_pres_split: the string, split into tokens around spaces,
                             but with the spaces conserved as tokens.
  """
  i = 0
  c = 0
  while (i < len(space_s_pres_split) and
         (space_s_pres_split[i].isspace() or
         (space_s_pres_split[i] == ""))):
    c += len(space_s_pres_split[i])
    i += 1
  return c


def _rm_get_names_from_header(parts):
  """
  get repeat and seq. name from repeatmasker alignment header line.

  An example header line is::

    239 29.42 1.92 0.97 chr1 11 17 (41) C XX#YY (74) 104 1 m_b1s502i1 4

  the genomic sequence name is always at position 4 (zero-based index); the
  name of the repeat is at position 9 if matching the reverse complement of
  the consensus sequence for the repeat and position 8 otherwise

  :param parts: the header line, as a tokenized list.
  :return: tuple of (name of genomic sequence, name of repeat sequence)
  """
  assert((parts[8] == "C" and len(parts) == 15) or (len(parts) == 14))
  return (parts[4], parts[8]) if len(parts) == 14 else (parts[4], parts[9])


def _rm_get_reference_coords_from_header(parts):
  """
  extract the reference (genomic sequence match) coordinates of a repeat
  occurrence from a repeatmakser header line. An example header line is::

    239 29.42 1.92 0.97 chr1 11 17 (41) C XX#YY (74) 104 1 m_b1s502i1 4

  the genomic start and end are always at positions 5 and 6 resepctively. In
  the repeatmasker format, the end is inclusive, but in pyokit end coordinates
  are exclusive, so we adjust it when we parse here.

  :param parts: the header line, as a tokenized list.
  :return: tuple of (start, end)
  """
  s = int(parts[5])
  e = int(parts[6]) + 1
  if (s >= e):
    raise AlignmentIteratorError("invalid repeatmakser header: " +
                                 " ".join(parts))
  return (s, e)


def _rm_get_repeat_coords_from_header(parts):
  """
  extract the repeat coordinates of a repeat masker match from a header line.

  An example header line is::

    239 29.42 1.92 0.97 chr1 11 17 (41) C XX#YY (74) 104 1 m_b1s502i1 4
    239 29.42 1.92 0.97 chr1 11 17 (41) XX#YY 1 104 (74) m_b1s502i1 4

  if the match is to the reverse complement, the start and end coordinates are
  at positions 11 and 12 (zero-based indexes), otherwise they're at positions
  9 and 10. In the later case, the 'start' is the earlier number and the end
  is the larger one. In reverse complement matches, RM lists the 'start' as the
  larger number and the end as the smaller one. We swap these around to match
  the Pyokit convention of start < end always and also adjust the end so it is
  not inclusive of the last position

  :param parts: the header line, as a tokenized list.
  :return: tuple of (start, end)
  """
  assert((parts[8] == "C" and len(parts) == 15) or (len(parts) == 14))
  if len(parts) == 14:
    s = int(parts[9])
    e = int(parts[10]) + 1
  else:
    s = int(parts[12])
    e = int(parts[11]) + 1
  if (s >= e):
    raise AlignmentIteratorError("invalid repeatmakser header: " +
                                 " ".join(parts))
  return (s, e)


def _rm_is_reverse_comp_match(parts):
  """
  determine whether a repeat occurrence is a match to the reverse complement
  of the concensus. Headers look like this::

    239 29.42 1.92 0.97 chr1 11 17 (41) C XX#YY (74) 104 1 m_b1s502i1 4

  If the match is to the reverse complement, then there is a "C" at position
  8 (zero-based index) and a total of 15 fields; otherwise the "C" is missing
  and there are only 14 fields.

  :param parts: the header line, as a tokenized list.
  """
  assert((parts[8] == "C" and len(parts) == 15) or (len(parts) == 14))
  return len(parts) == 15


def _rm_get_remaining_genomic_from_header(parts):
  """
  get the remaining number of bases that are on the genomic sequence after
  the match from the header line of a repeatmasker alignment. An example header
  line is::

    239 29.42 1.92 0.97 chr1 11 17 (41) C XX#YY (74) 104 1 m_b1s502i1 4

  The remaining genomic bases are always at position 7 (zero-based index)
  """
  return int(parts[7][1:-1])


def _rm_get_remaining_repeat_from_header(parts):
  """
  get the remaining number of bases that are on the repeat consensus after
  the match from the header line of a repeatmasker alignment. An example header
  line is::

    239 29.42 1.92 0.97 chr1 11 17 (41) XX#YY 1 104 (74) m_b1s502i1 4
    239 29.42 1.92 0.97 chr1 11 17 (41) C XX#YY (74) 104 1 m_b1s502i1 4

  If the match is to the consensus, this number is at position 11 (zero-based
  index), while a match to the reverse complement places it at position 10.
  The parenthese indicate it is a negative strand coordinate.
  """
  if _rm_is_reverse_comp_match(parts):
    return int(parts[10][1:-1])
  else:
    return int(parts[11][1:-1])


def _rm_parse_header_line(parts, meta_data):
  """
  parse a repeatmasker alignment header line and place the extracted meta-data
  into the provided dictionary. An example header line is::

    239 29.42 1.92 0.97 chr1 11 17 (41) C XX#YY (74) 104 1 m_b1s502i1 4

  If the alignment is to the consensus, this will have 14 fields; to the
  reverse complement of the repeat consensus and it'll have 15. Fields as
  follows:

  =====  ========    ==========================================================
  Num    Value       description
  =====  ========    ==========================================================
  0      239         Smith-Waterman score of the match
  1      29.42       percent substitutions in match compared to the consensus
  2      1.92        percent bases opposite gap in query seq (deleted bp)
  3      0.97        percent bases opposite gap in repeat consensus (insert bp)
  4      chr1        The name of the reference sequence (usually genomic chrom)
  5      11          Start location in reference; always pos strand, I think..
  6      17          End location in reference; always pos strand, I think..
  7      (41)        Distance to end of ref (alt., neg strand co-ord of end)
  8      C           Alignment is to rev. comp. of consensus. Ommited otherwise
  8/9    XX#YY       XX is the name of the repeat and YY is the family
  9/10   (74)        if +'ve, strt coord; else bases before start on neg strand
  10/11  104         if +'ve, end coord; else start of match in pos strand
  11/12  1           if +'ve, num bases aftr mtch to end; else end in + coords
  12/13  m_b1s5..    ....?
  13/14  4           unique ID
  =====  ========    ==========================================================

  Note that repeat-masker calls the larger coordinate the start when the match
  is to the reverse complement; we swap these internally so start < end always,
  regardless of whether the match is to the consensus or the reverse complement
  of the consensus

  Each field is mapped to a key as follows

  =====  =================================================================
  Num    Key (these are defined in pyokit.datastruct.multipleAlignment.py)
  =====  =================================================================
  0      ALIG_SCORE_KEY
  1      PCENT_SUBS_KEY
  2      PCENT_S1_INDELS_KEY
  3      PCENT_S2_INDELS_KEY
  4      S1_NAME_KEY
  5      S1_START_KEY
  6      S1_END_KEY
  7      S1_END_NEG_STRAND_KEY
  8      S2_REVERSE_COMP_KEY
  8/9    S2_NAME_KEY
  9/10   S2_START_KEY (S2 + strand) / S2_START_NEG_STRAND_KEY (S2 - strand)
  10/11  S2_END_KEY (S2 + strand) / S2_START_KEY (S2 - strand)
  11/12  S2_END_NEG_STRAND_KEY (S2 + strand) / S2_END_KEY (S2 - strand)
  12/13  UNKNOWN_RM_HEADER_FIELD_KEY
  13/14  RM_ID_KEY
  =====  =================================================================

  :param parts:     the header line, as a tokenized list.
  :param meta_data: dictionary; resultant key-value pairs placed into this.
  """
  meta_data[ALIG_SCORE_KEY] = parts[0]
  meta_data[PCENT_SUBS_KEY] = float(parts[1])
  meta_data[PCENT_S1_INDELS_KEY] = float(parts[2])
  meta_data[PCENT_S2_INDELS_KEY] = float(parts[3])
  meta_data[ANNOTATION_KEY] = ""

  if parts[8] == "C":
    meta_data[UNKNOWN_RM_HEADER_FIELD_KEY] = parts[13]
    meta_data[RM_ID_KEY] = int(parts[14])
  else:
    meta_data[UNKNOWN_RM_HEADER_FIELD_KEY] = parts[12]
    meta_data[RM_ID_KEY] = int(parts[13])


def _rm_name_match(s1, s2):
  """
  determine whether two sequence names from a repeatmasker alignment match.

  :return: True if they are the same string, or if one forms a substring of the
           other, else False
  """
  m_len = min(len(s1), len(s2))
  return s1[:m_len] == s2[:m_len]


def _rm_parse_meta_line(parts):
  p_locs = []
  for i in range(0, len(parts)):
    if parts[i].strip() == "=":
      p_locs.append(i)
  if len(p_locs) != 1:
    return ROUNDTRIP_KEY, " ".join(parts)
  else:
    k = " ".join(parts[:p_locs[0]])
    v = " ".join(parts[p_locs[0] + 1:])
    return k.strip(), v.strip()


def _rm_extract_sequence_and_name(alig_str_parts, s1_name, s2_name):
  """
  parse an alignment line from a repeatmasker alignment and return the name
  of the sequence it si from and the sequence portion contained in the line.

  :param alig_str_parts: the alignment string, split around whitespace as list
  :param s1_name: the name of the first sequence in the alignment this line is
                  from
  :param s2_name: the name of the second sequence in the alignment this line is
                  from
  :return: a tuple of name and sequence string; name will always be either
           s1_name or s2_name
  :raise AlignmentIteratorError: if the line doesn't have the expected number
                                 of elements, or the name does not match
                                 either of s1_name or s2_name
  """
  # first, based on the number of parts we have we'll guess whether its a
  # reverse complement or not
  if len(alig_str_parts) == 4:
    # expect the first element to amtch something..
    nm = alig_str_parts[0]
    seq = alig_str_parts[2]
  elif len(alig_str_parts) == 5:
    # expect the second element to match something...
    nm = alig_str_parts[1]
    seq = alig_str_parts[3]
  else:
    raise AlignmentIteratorError("failed parsing alignment line '" +
                                 " ".join(alig_str_parts) + "'; reason: " +
                                 "expected this line to have 4 or 5 " +
                                 "elements, but it has " +
                                 str(len(alig_str_parts)))
  if _rm_name_match(nm, s1_name):
    return s1_name, seq
  elif _rm_name_match(nm, s2_name):
    return s2_name, seq
  else:
    raise AlignmentIteratorError("failed parsing alignment line '" +
                                 " ".join(alig_str_parts) + "'; reason: " +
                                 "extracted alignment name (" + nm + ") " +
                                 "did not match either sequence name from " +
                                 "header line (" + s1_name + " or " +
                                 s2_name + ")")


def repeat_masker_alignment_iterator(fn, index_friendly=True, verbose=False):
  """
  Iterator for repeat masker alignment files; yields multiple alignment objects.

  Iterate over a file/stream of full repeat alignments in the repeatmasker
  format. Briefly, this format is as follows: each record (alignment) begins
  with a header line (see _rm_parse_header_line documentation for details of
  header format), followed by the alignment itself (example below) and finally
  a set of key-value meta-data pairs.

  The actual alignment looks like this::

    chr1               11 CCCTGGAGATTCTTATT--AGTGATTTGGGCT 41
                             ii        v   -- v  i i    v
    C MER5B#DNA/hAT    10 CCCCAGAGATTCTGATTTAATTGGTCTGGGGT 42

    chr1               42 GACTG 47
                           v
    C MER5B#DNA/hAT    43 CACTG 48

  The 'C' indicates that its the reverse complement of the consensus. The
  central string gives information about matches; "-" indicates an
  insertion/deletion, "i" a transition (G<->A, C<->T) and "v" a transversion
  (all other substitutions).

  :param fh:             filename or stream-like object to read from.
  :param index_friendly: if True, we will ensure the file/stream
                         position is before the start of the record when we
                         yield it; this requires the ability to seek within
                         the stream though, so if iterating over a
                         stream wtihout that ability, you'll have to set this
                         to false. Further, this will disable buffering for
                         the file, to ensure file.tell() behaves correctly,
                         so a performance hit will be incurred.
  :param verbose:        if true, output progress messages to stderr.
  """
  # step 1 -- build our iterator for the stream..
  try:
    fh = open(fn)
  except (TypeError):
    fh = fn
  iterable = fh
  if index_friendly:
    iterable = iter(fh.readline, '')

  # build progress indicator, if we want one and we're able to
  if verbose:
    try:
      m_fn = ": " + fh.name
    except TypeError:
      m_fn = ""
    try:
      current = fh.tell()
      fh.seek(0, 2)
      total_progress = fh.tell()
      fh.seek(current)
      pind = ProgressIndicator(totalToDo=total_progress,
                               messagePrefix="completed",
                               messageSuffix="of processing repeat-masker "
                                             "alignment file" + m_fn)
    except IOError:
      pind = None

  old_fh_pos = None
  new_fh_pos = fh.tell()

  s1 = None
  s2 = None
  s1_name = None
  s2_name = None
  s1_start = None
  s1_end = None
  s2_start = None
  s2_end = None
  meta_data = None
  alignment_line_counter = 0
  alig_l_space = 0
  prev_seq_len = 0
  rev_comp_match = None
  remaining_repeat = None
  remaining_genomic = None

  for line in iterable:
    if verbose and pind is not None:
      pind.done = fh.tell()
      pind.showProgress()

    if index_friendly:
      old_fh_pos = new_fh_pos
      new_fh_pos = fh.tell()
    line = line.rstrip()
    if line.lstrip() == "" and alignment_line_counter % 3 != 1:
      continue

    s_pres_split = re.split(r'(\s+)', line)
    parts = [x for x in s_pres_split if not (x.isspace() or x == "")]

    n = len(parts)
    for i in REPEATMASKER_FIELDS_TO_TRIM:
      if n >= i + 1:
        parts[i] = parts[i].strip()

    # decide what to do with this line -- is it a header line, part of the
    # alignment or a meta-data key-value line
    if alignment_line_counter % 3 == 1:
      if (REPEATMASKER_VALIDATE_MUTATIONS and
         not _rm_is_valid_annotation_line(line)):
        raise IOError("invalid mutation line: " + line)
      l_space = _rm_compute_leading_space(s_pres_split) - alig_l_space
      pad_right = prev_seq_len - (l_space + len(line.strip()))
      meta_data[ANNOTATION_KEY] += ((' ' * l_space) + line.strip() +
                                    (' ' * pad_right))
      alignment_line_counter += 1
    elif _rm_is_header_line(parts, n):
      if not (s1 is None and s2 is None and meta_data is None):
        if ANNOTATION_KEY in meta_data:
          meta_data[ANNOTATION_KEY] = meta_data[ANNOTATION_KEY].rstrip()
        if index_friendly:
          fh.seek(old_fh_pos)
        ss1 = Sequence(s1_name, s1, s1_start, s1_end, "+", remaining_genomic)
        s2s = "-" if rev_comp_match else "+"
        ss2 = Sequence(s2_name, s2, s2_start, s2_end, s2s, remaining_repeat)
        yield PairwiseAlignment(ss1, ss2, meta_data)
        if index_friendly:
          fh.seek(new_fh_pos)
      meta_data = {}
      s1 = ""
      s2 = ""
      s1_name, s2_name = _rm_get_names_from_header(parts)
      s1_start, s1_end = _rm_get_reference_coords_from_header(parts)
      s2_start, s2_end = _rm_get_repeat_coords_from_header(parts)
      rev_comp_match = _rm_is_reverse_comp_match(parts)
      remaining_repeat = _rm_get_remaining_repeat_from_header(parts)
      remaining_genomic = _rm_get_remaining_genomic_from_header(parts)

      _rm_parse_header_line(parts, meta_data)
      alignment_line_counter = 0
    elif _rm_is_alignment_line(parts, s1_name, s2_name):
      alignment_line_counter += 1
      name, seq = _rm_extract_sequence_and_name(parts, s1_name, s2_name)
      if name == s1_name:
        s1 += seq
      elif name == s2_name:
        s2 += seq
      alig_l_space = _rm_compute_leading_space_alig(s_pres_split, seq)
      prev_seq_len = len(seq)
    else:
      k, v = _rm_parse_meta_line(parts)
      meta_data[k] = v
  if index_friendly:
    fh.seek(old_fh_pos)
  ss1 = Sequence(s1_name, s1, s1_start, s1_end, "+", remaining_genomic)
  s2s = "-" if rev_comp_match else "+"
  ss2 = Sequence(s2_name, s2, s2_start, s2_end, s2s, remaining_repeat)
  yield PairwiseAlignment(ss1, ss2, meta_data)
  if index_friendly:
    fh.seek(new_fh_pos)


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestAlignmentIterators(unittest.TestCase):

  def setUp(self):
    # set up a repeat-masker file...
    alig_1_header = "283 26.37 4.21 0.00 chr1 15 67 (266) C A#B (119) " +\
                    "141 85 m_b1s601i0 5                              "
    alig_1 = "  chr1  15 CCACTGTACA-ATGGGGAAACT--GGCCC 40     \n" +\
             "              i v    -   i       --   v         \n" +\
             "C A#B  141 CCATTTTACAGATGAGGAAACTGAGGCAC 113    \n" +\
             "                                                \n" +\
             "  chr1  41 AGAGCAAGGCAAAAGCAGCGCTGGG-TA 67      \n" +\
             "           v   v  vv ivi    v  i    - v         \n" +\
             "C A#B  112 CAGCTAGTAAGTGGCAGAGCCGGGATTC 85        "
    alig_1_m = "Matrix = 25p47g.matrix                       \n" +\
               "Kimura (with divCpGMod) = 29.95              \n" +\
               "Transitions / transversions = 1.40 (14/10)   \n" +\
               "Gap_init rate = 0.03 (3 / 90), avg. gap size = 1.33 (4 / 3)"

    alig_2_header = "318 22.50 3.61 0.00 chr1 15266 15323 (249235276) C " +\
                    "MIR3#SINE/MIR (65) 143 84 m_b1s601i1 10"
    alig_2 = "  chr1          15266 GAAACT--GGCCCAGAGAGGTGAGGCAGCG 15293 \n" +\
             "                            --               i iii         \n" +\
             "C MIR3#SINE/MIR   143 GAAACTGAGGCCCAGAGAGGTGAAGTGACG 114   \n" +\
             "                                                           \n" +\
             "  chr1          15294 GGTCACAGAGCAAGGCAAAAGCGCGCTGGG 15323 \n" +\
             "                             v   ?  vi ivi    v            \n" +\
             "C MIR3#SINE/MIR   113 GGTCACACAGCKAGTTAGTGGCGAGCTGGG 84"
    alig_2_m = "Matrix = 25p47g.matrix                                   \n" +\
               "Kimura (with divCpGMod) = 26.25                          \n" +\
               "Transitions / transversions = 2.40 (12/5)                \n" +\
               "Gap_init rate = 0.03 (2 / 79), avg. gap size = 1.50 (3 / 2)"

    alig_3_header = "18 23.18 0.00 1.96 chr1 15798 15830 (249234772) " +\
                    "(TGCTCC)n#Simple_repeat 1 32 (0) m_b1s252i0 15"
    alig_3 = "  chr1          15798 GCTGCTTCTCCAGCTTTCGCTCCTTCATGCT 15828 \n" +\
             "                         v  v    v   iii      v - v         \n" +\
             "  (TGCTCC)n#Sim     1 GCTCCTGCTCCTGCTCCTGCTCCTGC-TCCT 30    \n" +\
             "                                                            \n" +\
             "  chr1          15829 GC 15830                              \n" +\
             "                                                            \n" +\
             "  (TGCTCC)n#Sim    31 GC 32                                    "
    alig_3_m = "Matrix = Unknown                                   \n" +\
               "Transitions / transversions = 0.43 (3/7)           \n" +\
               "Gap_init rate = 0.02 (1 / 51), avg. gap size = 1.00 (1 / 1)"

    alig_4_header = "487 20.75 0.93 0.93 chr1 158389 158409 (249092126) C " +\
                    "Charlie29b#DNA/hAT-Charlie (532) 662 641 m_b3s502i21 231"
    alig_4 = "  chr1          158389 TAGAATTTTTGTGGCAT-ATGA 158409    \n" +\
             "                          i ii v ii     -  vi           \n" +\
             "C Charlie29b#DN    662 TAAAGCTGGGCGTTATTGATGA 641       \n"
    alig_4_m = ""

    self.rm_tc1_records = [alig_1_header + "\n\n" + alig_1 + "\n\n" + alig_1_m,
                           alig_2_header + "\n\n" + alig_2 + "\n\n" + alig_2_m,
                           alig_3_header + "\n\n" + alig_3 + "\n\n" + alig_3_m,
                           alig_4_header + "\n\n" + alig_4 + "\n\n" + alig_4_m]
    self.rm_rc_1_input = "\n\n".join(self.rm_tc1_records)

  def test_rm_iter_one_part_ann_line(self):
    """
    This test for the repeatmakser iterator has an annotation line with just a
    single element in it; this exposes a bug that previously existed (now
    fixed) where the lines are expected to have >= 2 elements.
    """
    debug = False

    # set up test
    alig_header = "10304 12.32 4.41 4.46 chrUn_gl000247 36396 36422 (0) " +\
                  "MER57-int#LTR/ERV1 5219 5245 (1398) m_b1s701i6 4562661"
    alig = "  chrUn_gl00024 36396 TTAATGTGAACAGCTTTTCCCAAGATC 36422\n" +\
           "                        i                              \n" +\
           "  MER57-int#LTR  5219 TTGATGTGAACAGCTTTTCCCAAGATC 5245"
    test_in = alig_header + "\n\n" + alig

    # run test and collect results
    results = [r for r in
               repeat_masker_alignment_iterator(StringIO.StringIO(test_in))]
    self.failUnlessEqual(len(results), 1)
    rm_str = _to_repeatmasker_string(results[0], m_name_width=13)

    # check results
    alig_actual = [x for x in map(str.rstrip, test_in.split("\n"))
                   if x.strip() != ""]
    alig_result = [x for x in map(str.rstrip, rm_str.split("\n"))
                   if x.strip() != ""]
    if debug:
      print ""
      print "expected: " + str(alig_actual)
      print "got:      " + str(alig_result)
    self.failUnlessEqual(alig_actual, alig_result)

  def test_repeat_masker_alignment_iterator(self):
    """Test roundtrip of repeatmasker alignment."""
    debug = False
    s_io = StringIO.StringIO(self.rm_rc_1_input)
    alig_iter = repeat_masker_alignment_iterator(s_io)
    results = [r for r in alig_iter]
    self.failUnlessEqual(len(results), len(self.rm_tc1_records))
    for i, trail_meta_size, c_width, m_width in [(0, 4, 29, None),
                                                 (1, 4, 30, None),
                                                 (2, 3, 31, 13),
                                                 (3, 0, 22, 13)]:
      rm_str = _to_repeatmasker_string(results[i], column_width=c_width,
                                       m_name_width=m_width)
      if debug:
        sys.stderr.write("===============================\n")
        sys.stderr.write(rm_str + "\n")
        sys.stderr.write("*******************************\n")
        sys.stderr.write(self.rm_tc1_records[i] + "\n")
        sys.stderr.write("===============================\n")

      # strip out the last few lines; these should all be their, but order
      # isn't important.
      alig_actual = [x for x in map(str.rstrip,
                     self.rm_tc1_records[i].split("\n")[:-trail_meta_size])
                     if x.strip() != ""]
      meta_actual = map(str.rstrip,
                        self.rm_tc1_records[i].split("\n")[-trail_meta_size:])
      alig_result = [x for x in map(str.rstrip,
                                    rm_str.split("\n")[:-trail_meta_size])
                     if x.strip() != ""]
      meta_result = map(str.rstrip, rm_str.split("\n")[-trail_meta_size:])

      if debug:
        sys.stderr.write(str(alig_actual) + "\n")
        sys.stderr.write(str(alig_result) + "\n")
      self.failUnlessEqual(alig_actual, alig_result)
      self.failUnlessEqual(set(meta_actual), set(meta_result))

  def test_repeat_masker_on_demand_load(self):
    """
    Tests wrapping the alignment iterator in an index and using this index
    to build RM alignment objects that are loaded on-demand from the indexed
    stream.
    """
    from pyokit.io.indexedFile import IndexedFile

    def extract_UID(rm_alignment):
      return rm_alignment.meta[RM_ID_KEY]

    s_io = StringIO.StringIO(self.rm_rc_1_input)
    index = IndexedFile(s_io, repeat_masker_alignment_iterator, extract_UID)

    for i, trail_meta_size, c_width, m_width, rm_id in [(0, 4, 29, None, 5),
                                                        (1, 4, 30, None, 10),
                                                        (2, 3, 31, 13, 15),
                                                        (3, 0, 22, 13, 231)]:
      on_d_alig = JustInTimePairwiseAlignment(index, rm_id)
      on_d_str = _to_repeatmasker_string(on_d_alig, column_width=c_width,
                                         m_name_width=m_width)

      # strip out the last few lines; these should all be their, but order
      # isn't important.
      alig_actual = [x for x in map(str.rstrip,
                     self.rm_tc1_records[i].split("\n")[:-trail_meta_size])
                     if x.strip() != ""]
      meta_actual = map(str.rstrip,
                        self.rm_tc1_records[i].split("\n")[-trail_meta_size:])
      alig_result = [x for x in map(str.rstrip,
                                    on_d_str.split("\n")[:-trail_meta_size])
                     if x.strip() != ""]
      meta_result = map(str.rstrip, on_d_str.split("\n")[-trail_meta_size:])

      self.failUnlessEqual(alig_actual, alig_result)
      self.failUnlessEqual(set(meta_actual), set(meta_result))


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == '__main__':
    unittest.main()
