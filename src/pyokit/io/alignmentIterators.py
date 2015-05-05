"""
  Date of Creation: 11th Dec 2014

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
import sys
import StringIO
import unittest
import re

# Pyokit imports
from pyokit.util.progressIndicator import ProgressIndicator
from pyokit.datastruct.multipleAlignment import PairwiseAlignment
from pyokit.datastruct import multipleAlignment
from pyokit.datastruct.multipleAlignment import JustInTimePairwiseAlignment


###############################################################################
#                                 CONSTANTS                                   #
###############################################################################

REPEATMASKER_FIELDS_TO_TRIM = [0, 1, 4, 7, 8, 9, 11]
REPEATMASKER_VALIDATE_MUTATIONS = True
REPEATMASKER_VALID_ANN_CHARS = ['-', 'i', 'v', ' ', '?']


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
  returns true if the tokenized line is a repeatmasker alignment line

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
  headers have no special structure or symbol to mark them...
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
  meta_data[multipleAlignment.ALIG_SCORE_KEY] = parts[0]
  meta_data[multipleAlignment.PCENT_SUBS_KEY] = float(parts[1])
  meta_data[multipleAlignment.PCENT_S1_INDELS_KEY] = float(parts[2])
  meta_data[multipleAlignment.PCENT_S2_INDELS_KEY] = float(parts[3])
  meta_data[multipleAlignment.ANNOTATION_KEY] = ""
  meta_data[multipleAlignment.S1_NAME_KEY] = parts[4]
  meta_data[multipleAlignment.S1_START_KEY] = int(parts[5])
  meta_data[multipleAlignment.S1_END_KEY] = int(parts[6]) + 1
  meta_data[multipleAlignment.S1_END_NEG_STRAND_KEY] = int(parts[7][1:-1])

  if parts[8] == "C" :
    meta_data[multipleAlignment.S2_REVERSE_COMP_KEY] = True
    meta_data[multipleAlignment.S2_NAME_KEY] = parts[9]
    meta_data[multipleAlignment.S2_START_NEG_STRAND_KEY] = int(parts[10][1:-1])
    meta_data[multipleAlignment.S2_START_KEY] = int(parts[12])
    meta_data[multipleAlignment.S2_END_KEY] = int(parts[11]) + 1
    meta_data[multipleAlignment.UNKNOWN_RM_HEADER_FIELD_KEY] = parts[13]
    meta_data[multipleAlignment.RM_ID_KEY] = int(parts[14])
  else:
    meta_data[multipleAlignment.S2_NAME_KEY] = parts[8]
    meta_data[multipleAlignment.S2_START_KEY] = int(parts[9])
    meta_data[multipleAlignment.S2_END_KEY] = int(parts[10]) + 1
    meta_data[multipleAlignment.S2_END_NEG_STRAND_KEY] = int(parts[11][1:-1])
    meta_data[multipleAlignment.UNKNOWN_RM_HEADER_FIELD_KEY] = parts[12]
    meta_data[multipleAlignment.RM_ID_KEY] = int(parts[13])

  # internally, we always require start < end
  assert(meta_data[multipleAlignment.S1_START_KEY] <
         meta_data[multipleAlignment.S1_END_KEY])
  assert(meta_data[multipleAlignment.S2_START_KEY] <
         meta_data[multipleAlignment.S2_END_KEY])


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
    if parts[i].strip() == "=" :
      p_locs.append(i)
  if len(p_locs) != 1:
    return multipleAlignment.ROUNDTRIP_KEY, " ".join(parts)
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
    raise AlignmentIteratorError("failed parsing alignment line '"
                                 + " ".join(alig_str_parts) + "'; reason: "
                                 + "expected this line to have 4 or 5 "
                                 + "elements, but it has "
                                 + str(len(alig_str_parts)))
  if _rm_name_match(nm, s1_name):
    return s1_name, seq
  elif _rm_name_match(nm, s2_name):
    return s2_name, seq
  else:
    raise AlignmentIteratorError("failed parsing alignment line '"
                                 + " ".join(alig_str_parts) + "'; reason: "
                                 + "extracted alignment name (" + nm + ") "
                                 + "did not match either sequence name from "
                                 + "header line (" + s1_name + " or "
                                 + s2_name + ")")


def repeat_masker_alignment_iterator(fn, index_friendly=True, verbose=False):
  """
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
    try :
      m_fn = ": " + fh.name
    except TypeError:
      m_fn = ""
    try :
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
  meta_data = None
  alignment_line_counter = 0
  alig_l_space = 0
  prev_seq_len = 0

  for line in iterable:
    if verbose and pind != None :
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
      meta_data[multipleAlignment.ANNOTATION_KEY] += ((' ' * l_space)
                                                      + line.strip()
                                                      + (' ' * pad_right))
      alignment_line_counter += 1
    elif _rm_is_header_line(parts, n):
      if not (s1 == None and s2 == None and meta_data == None):
        if multipleAlignment.ANNOTATION_KEY in meta_data:
          meta_data[multipleAlignment.ANNOTATION_KEY] = \
              meta_data[multipleAlignment.ANNOTATION_KEY].rstrip()
        if index_friendly:
          fh.seek(old_fh_pos)
        yield PairwiseAlignment(s1, s2, meta_data)
        if index_friendly:
          fh.seek(new_fh_pos)
      meta_data = {}
      s1 = ""
      s2 = ""
      _rm_parse_header_line(parts, meta_data)
      alignment_line_counter = 0
    elif _rm_is_alignment_line(parts, meta_data[multipleAlignment.S1_NAME_KEY],
                               meta_data[multipleAlignment.S2_NAME_KEY]):
      alignment_line_counter += 1
      s1_n = meta_data[multipleAlignment.S1_NAME_KEY]
      s2_n = meta_data[multipleAlignment.S2_NAME_KEY]
      name, seq = _rm_extract_sequence_and_name(parts, s1_n, s2_n)
      if name == s1_n:
        s1 += seq
      elif name == s2_n:
        s2 += seq
      alig_l_space = _rm_compute_leading_space_alig(s_pres_split, seq)
      prev_seq_len = len(seq)
    else:
      k, v = _rm_parse_meta_line(parts)
      meta_data[k] = v
  if index_friendly:
    fh.seek(old_fh_pos)
  yield PairwiseAlignment(s1, s2, meta_data)
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
    alig_3 = "  chr1          15798 GCTGCTTCTCCAGCTTTCGCTCCTTCATGCT 15828  \n" +\
             "                         v  v    v   iii      v - v          \n" +\
             "  (TGCTCC)n#Sim     1 GCTCCTGCTCCTGCTCCTGCTCCTGC-TCCT 30     \n" +\
             "                                                             \n" +\
             "  chr1          15829 GC 15830                               \n" +\
             "                                                             \n" +\
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
    rm_str = results[0].to_repeat_masker_string(m_name_width=13)

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
    """
    This is a roundtrip test.
    """

    debug = False
    s_io = StringIO.StringIO(self.rm_rc_1_input)
    alig_iter = repeat_masker_alignment_iterator(s_io)
    results = [r for r in alig_iter]
    self.failUnlessEqual(len(results), len(self.rm_tc1_records))
    for i, trail_meta_size, c_width, m_width in [(0, 4, 29, None),
                                                 (1, 4, 30, None),
                                                 (2, 3, 31, 13),
                                                 (3, 0, 22, 13)]:
      rm_str = results[i].to_repeat_masker_string(column_width=c_width,
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
      return rm_alignment.meta[multipleAlignment.RM_ID_KEY]

    s_io = StringIO.StringIO(self.rm_rc_1_input)
    index = IndexedFile(s_io, repeat_masker_alignment_iterator, extract_UID)

    for i, trail_meta_size, c_width, m_width, rm_id in [(0, 4, 29, None, 5),
                                                        (1, 4, 30, None, 10),
                                                        (2, 3, 31, 13, 15),
                                                        (3, 0, 22, 13, 231)]:
      on_d_alig = JustInTimePairwiseAlignment(index, rm_id)
      on_d_str = on_d_alig.to_repeat_masker_string(column_width=c_width,
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
