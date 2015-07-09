"""
Date of Creation: 11th Dec 2014.

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

# backported enum package
from enum import Enum

# pyokit imports
from pyokit.datastruct.sequence import Sequence
from pyokit.datastruct.sequence import UnknownSequence
from pyokit.datastruct.sequence import InvalidSequenceCoordinatesError
from pyokit.datastruct.sequence import GAP_CHAR
from pyokit.util.meta import decorate_all_methods_and_properties
from pyokit.util.meta import just_in_time_method, just_in_time_property
from pyokit.common.pyokitError import PyokitError


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class MultipleAlignmentError(PyokitError):

  """Errors that occur when manipulating pairwise or multiple alig. objects."""

  def __init__(self, msg):
    """Constructor for MSA errors."""
    self.value = msg

  def __str__(self):
    """:return: string representation of this exception."""
    return repr(self.value)


class InvalidSequenceError(MultipleAlignmentError):

  """Throw when provided an invalid sequence index into the MSA."""

  def __init__(self, msg):
    """Constructor for InvalidSequenceErrors."""
    self.value = msg

  def __str__(self):
    """:return: string representation of this exception."""
    return repr(self.value)


class InvalidAlignmentCoordinatesError(MultipleAlignmentError):

  """Thrown when providing invalid coords. into an MSA or a sequence in it."""

  def __init__(self, msg):
    """Constructor for InvalidAlignmentCoordinatesErrors."""
    self.value = msg

  def __str__(self):
    """:return: string representation of this exception."""
    return repr(self.value)


###############################################################################
#                                ENUMERATIONS                                 #
###############################################################################
class MissingSequenceHandler(Enum):

  """An enum of ways to handle missing sequences."""

  SKIP = 1
  TREAT_AS_ALL_GAPS = 2


###############################################################################
#                          MULTIPLE ALIGNMENT CLASS                           #
###############################################################################

class MultipleSequenceAlignment(object):

  """
  An alignment of two or more sequences.

  :param sequences: a list of sequence objects; all have to be the same length
  """

  def __init__(self, sequences, meta_data=None, permit_one_seq=True):
    """Constructor for MSAs; see class docstring for param descriptions."""
    # must get at least 2 sequences
    if len(sequences) == 0:
      msg = "Constructing multiple alignment object failed; no sequences"
      raise MultipleAlignmentError(msg)
    if len(sequences) == 1 and not permit_one_seq:
      msg = "Constructing multiple alignment object failed; expected at " +\
            "least 2 sequences, found one sequence: " +\
            str(sequences[0].name) + " --> " + str(sequences[0].sequenceData)
      raise MultipleAlignmentError(msg)

    # make sure they're all the same length
    lengths = set()
    for s in sequences:
      if isinstance(s, UnknownSequence):
        continue
      lengths.add(len(s))
    if len(lengths) != 1:
      msg = "Invalid multiple alignment, sequences have different lengths: " +\
            ", ".join([str(x) for x in lengths]) + "; sequences are " +\
            ",".join(str(x.sequenceData) for x in sequences)
      raise MultipleAlignmentError(msg)
    self.length = list(lengths)[0]

    # internally we store the sequences in a dictionary indexed by seq name
    self.sequences = {}
    for s in sequences:
      self.sequences[s.name] = s

    self._meta = meta_data

  def __getitem__(self, seq_name):
    """:return: the sequence with the given name in this alignment."""
    try:
      return self.sequences[seq_name]
    except KeyError:
      raise InvalidSequenceError("No such sequence in alignemnt: " +
                                 str(seq_name) + "; sequences in alignment: " +
                                 ", ".join(self.sequences.keys()))

  def get_column(self, position, missing_seqs=MissingSequenceHandler.SKIP):
    """
    return a column from an alignment as a dictionary indexed by seq. name.

    :param position:     the index to extract; these are in alignment
                         co-ordinates, which are one-based, so the first column
                         has index 1, and the final column has
                         index == size(self).
    :param missing_seqs: how to treat sequence with no actual sequence data for
                         the column.
    :return: dictionary where keys are sequence names and values are
             nucleotides (raw strings).
    """
    res = {}
    for k in self.sequences:
      if isinstance(self.sequences[k], UnknownSequence):
        if missing_seqs is MissingSequenceHandler.TREAT_AS_ALL_GAPS:
          res[k] = "-"
        elif missing_seqs is MissingSequenceHandler.SKIP:
          continue
      else:
        res[k] = self.sequences[k][position - 1]
    return res

  def __iter__(self):
    """:return: an iterator for the alignment, which yields sequences."""
    return self.sequences.__iter__()

  @property
  def meta(self):
    """:return: dictionary of key-value pairs; the alignment's meta data."""
    return self._meta

  def __str__(self):
    """:return: a string representation of the alignment."""
    res = ""
    keys = self.sequences.keys()
    for i in range(0, len(keys)):
      if i != 0:
        res += "\n"
      k = keys[i]
      res += str(self.sequences[k])
    return res

  def size(self):
    """:return: the length (number of columns) in this multiple alignment."""
    return self.length

  def num_seqs(self):
    """:return: the length (number of sequences) in this alignment."""
    return len(self.sequences)

  def alignment_to_sequence_coords(self, seq_name, start, end, trim=False):
    """
    convert an interval in the alignmnet into co-ordinates in one of the
    sequences. Alignment intervals are inclusive of start, but not end. They
    are one-based. Hence the full alignment has coords [1, N+1), where N is the
    length of the alignment (number of columns). Sequence coords follow the
    same conventions: one-based, inclusive of start but not end.

    :param seq_name: which sequence are the start and end coords for?
    :param start:    start of the interval in alignment co-ordinates
    :param end:      end of the interval in alignment co-ordinates
    :param trim:     if true, trim coordinates that fall partially outside
                     the alignment
    :return: a tuple with the start and end sequence coords, in that order;
             None if the interval in the alignment defined by (start, end)
             contains only gaps in the specified sequence. start < end always,
             even if the sequence is reverse complement.
    :raises AlignmentError: if the sequence number specifies a sequence not
                            in this alignment, or the coordinates fall
                            entirely outside the alignment (or partially
                            outside and trim == false), or they are not a valid
                            interval (start >= end)
    """
    start = 1 if start < 1 and trim else start
    end = self.size() + 1 if (end > self.size() + 1) and trim else end

    if (start < 1 or end > (self.size() + 1)):
      msg = "Coordinates fall partially outside alignemnt: " + str(start) +\
            ", " + str(end)
      raise InvalidAlignmentCoordinatesError(msg)
    if (end < 1 or start > self.size() + 1):
      msg = "Coordinates fall entirely outside alignment: " + str(start) +\
            ", " + str(end)
      raise InvalidAlignmentCoordinatesError(msg)
    if (end <= start):
      msg = "Invalid alignment coordinates: " + str(start) + ", " + str(end)
      raise InvalidAlignmentCoordinatesError(msg)

    seq = self[seq_name]
    pos_strand = seq.is_positive_strand()
    non_gaps = 0
    r_start = None
    r_end = None
    l_start = 0 if pos_strand else self.size() - 1
    l_end = end - 1 if pos_strand else start - 2
    l_step = 1 if pos_strand else -1
    for i in range(l_start, l_end, l_step):
      if seq[i] != GAP_CHAR:
        non_gaps += 1
        if ((pos_strand and r_start is None and (i + 1) >= start) or
           (not pos_strand and r_start is None and (i + 1) < end)):
          r_start = non_gaps + seq.start - 1
    if r_start is None:
      # we went through the whole region and didn't find a single non-gap char
      return None
    r_end = non_gaps + seq.start
    return (r_start, r_end)

  def sequence_to_alignment_coords(self, seq_name, start, end, trim=False):
    """
    convert an interval in one of the sequences into an interval in the
    alignment. Alignment intervals are inclusive of start, but not end. They
    are one-based. Hence the full alignment has coords [1, N+1), where N is the
    length of the alignment (number of columns). Sequence coords follow the
    same conventions: one-based, inclusive of start but not end.

    :param seq_name: which sequence are the start and end coords for?
    :param start:    start of the interval in sequence co-ordinates
    :param end:      end of the interval in sequence co-ordinates
    :param trim:     if true, trim coordinates that fall partially outside
                     the sequence
    :raises AlignmentError: if coordinates fall entirely outside the
                            sequence, or partially outside and trim == false
    """
    # check for valid order of start/end
    if end <= start:
      raise InvalidSequenceCoordinatesError("invalid region: " + str(start) +
                                            ", " + str(end))
    seq = self[seq_name]
    s_start = seq.start
    s_end = seq.end
    pos_strand = seq.is_positive_strand()

    # check that the start and end coords are at least partially in the seq
    if start > s_end or end < s_start:
      msg = "Cannot convert " + str(start) + ", " + str(end) + " to " +\
            "alignment coordinates; falls fully outside of sequence " +\
            str(s_start) + ", " + str(s_end)
      raise InvalidSequenceCoordinatesError(msg)

    # trim overlap if that option is sepcified, otherwise complain if outside
    if trim:
      start = s_start if start < s_start else start
      end = s_end if end > s_end else end
    elif start < s_start or end > s_end:
      msg = "Cannot convert " + str(start) + ", " + str(end) + " to " +\
            "alignment coordinates; falls artially outside of seq. " +\
            str(s_start) + ", " + str(s_end)
      raise InvalidSequenceCoordinatesError(msg)

    num_gaps = 0
    num_non_gaps = 0
    res = []
    current_start = None
    current_end = None
    l_start = 0 if pos_strand else self.size() - 1
    l_end = self.size() if pos_strand else 0
    l_step = 1 if pos_strand else -1
    for i in range(l_start, l_end, l_step):
      if seq[i] == GAP_CHAR:
        num_gaps += 1
      else:
        num_non_gaps += 1

      if num_non_gaps > end - s_start:    # done, past the end of the ROI
        break
      if num_non_gaps > start - s_start:  # within ROI still
        if seq[i] != GAP_CHAR:
          if current_start is None and current_end is None:
            current_start = i
            current_end = i + 1
          else:
            if ((pos_strand and seq[i - 1] == GAP_CHAR) or
                (not pos_strand and seq[i + 1] == GAP_CHAR)):
              # is the start of a new non-gapped region...
              res.append((current_start + 1, current_end + 1))
              current_start = i
              current_end = i + 1
            if pos_strand and seq[i - 1] != GAP_CHAR:
              # is continuation of non-gapped region
              current_end += 1
            if not pos_strand and seq[i + 1] != GAP_CHAR:
              # is continuation of non-gapped region
              current_start -= 1
    res.append((current_start + 1, current_end + 1))
    return res

  def liftover(self, origin, dest, o_start, o_end, trim=False):
    """liftover interval in one seq. of this pairwise alignment to the other.

    :param origin:  name of the origin seq (seq the input coordinates are for)
    :param dest:    name of the dest. seq (seq the result will be for)
    :param o_start: start of the interval (in sequence co-ordinates) to lift.
    :param o_end:   end of the interval (in seq. coords) to lift.
    """
    alig_cols = self.sequence_to_alignment_coords(origin, o_start,
                                                  o_end, trim=trim)
    res = []
    for s, e in alig_cols:
      t = self.alignment_to_sequence_coords(dest, s, e)
      if t is None:
        continue
      res.append(t)
    return res


###############################################################################
#                          PAIRWISE ALIGNMENT CLASS                           #
###############################################################################

class PairwiseAlignment(MultipleSequenceAlignment):

  """An alignment of two sequences (DNA, RNA, protein...).

  :param s1: the first sequence, with gaps
  :param s2: the second sequence, with gaps
  :param meta_data: a dictionary with key-value pairs representing meta-data
                    about this alignment.
  """

  def __init__(self, s1, s2, meta_data=None):
    """Constructor for pairwise alig. see class dostring for param details."""
    MultipleSequenceAlignment.__init__(self, [s1, s2], meta_data)
    self.s1_name = s1.name
    self.s2_name = s2.name

  @property
  def s1(self):
    """:return the first sequences in the alignment."""
    return self[self.s1_name]

  @property
  def s2(self):
    """:return the second sequences in the alignment."""
    return self[self.s2_name]

  def __str__(self):
    """:return: a string representation of this pairwise alignment."""
    return self.to_repeat_masker_string()


###############################################################################
#         ON-DEMAND LOADING OF PAIRWISE ALIGNMENTS FROM INDEXED FILE          #
###############################################################################

@decorate_all_methods_and_properties(just_in_time_method,
                                     just_in_time_property)
class JustInTimePairwiseAlignment(PairwiseAlignment):

  """A pairwise alignment that is loaded just-in-time from some factory object.

  A common pattern would be to provide an IndexedFile as the factory object

  :param factory: any object that implements the subscript operator such that
                  it accept the key as a unique identifier and returns the
                  alignment
  :param key:     any object which uniquely identifies this pariwise alignment
                  to the factory; i.e. a hash key
  """

  def __init__(self, factory, key):
    """Constructor for JIT pairwise aligs; see class doc. for param details."""
    self.factory = factory
    self.key = key
    self.item = None


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################
class TestAlignments(unittest.TestCase):

  """Unit tests for the multiple alignment module."""

  def setUp(self):
    """Set up a few alignments to use in the tests."""
    s1 = Sequence("s1", "-TCGCGTAGC---CGC-TAGCTGATGCGAT-CTGA", 100, 129)
    s2 = Sequence("s2", "ATCGCGTAGCTAGCGCG-AGCTG---CGATGCT--", 1000, 1029)
    s3 = Sequence("s3", "ATCGCGTAGCTAGCGCG-AGCTG---CGATGCT--", 969, 998, "-")

    self.pa1 = PairwiseAlignment(s1, s2)
    self.pa2 = PairwiseAlignment(s1, s3)
    self.msa1 = MultipleSequenceAlignment([s1, s2, s3])

  def test_ungapped_length(self):
    """Test computing the ungapped length of sequences in the alignment."""
    self.assertEqual(self.pa1["s1"].ungapped_len,
                     len("TCGCGTAGCCGCTAGCTGATGCGATCTGA"))
    self.assertEqual(self.pa1["s1"].ungapped_len,
                     len("ATCGCGTAGCTAGCGCGAGCTGCGATGCT"))
    self.assertEqual(self.pa2["s1"].ungapped_len,
                     len("TCGCGTAGCCGCTAGCTGATGCGATCTGA"))
    self.assertEqual(self.pa2["s3"].ungapped_len,
                     len("ATCGCGTAGCTAGCGCGAGCTGCGATGCT"))

  def test_get_column(self):
    """Test extracting single columns from multiple and pairwise aligs."""
    self.assertEqual(self.msa1.get_column(1), {"s1": "-", "s2": "A",
                                               "s3": "A"})
    self.assertEqual(self.msa1.get_column(5), {"s1": "C", "s2": "C",
                                               "s3": "C"})
    self.assertEqual(self.pa1.get_column(11), {"s1": "-", "s2": "T"})

  def test_sequence_to_alig_coord(self):
    """Test converting sequence coords into alig. coords."""
    # CGTAGC---CGC
    # CGTAGCTAGCGC
    self.assertEqual(self.pa1.sequence_to_alignment_coords("s1", 103, 112),
                     [(5, 11), (14, 17)])
    # 977 --> G---CGAT <-- 972
    self.assertEqual(self.pa2.sequence_to_alignment_coords("s3", 972, 977),
                     [(27, 31), (23, 24)])
    # should throw an exception if we give coordinates that extend outside the
    # sequence...
    self.assertRaises(InvalidSequenceCoordinatesError,
                      self.pa1.sequence_to_alignment_coords, "s1", 95, 112)
    # .. but not if we ask them to be trimmed ...
    r = self.pa1.sequence_to_alignment_coords("s1", 95, 112, trim=True)
    self.assertEqual(r, [(2, 11), (14, 17)])
    # .. but still if they fall entirely outside the sequence ..
    self.assertRaises(InvalidSequenceCoordinatesError,
                      self.pa1.sequence_to_alignment_coords, "s1", 95, 99)
    # test lifting over coordinates that overlaps the end of the alignement
    r = self.pa1.sequence_to_alignment_coords("s1", 124, 130, trim=True)
    self.assertEqual(r, [(30, 31), (32, 36)])

  def test_alig_to_sequence_coords(self):
    """Test convert alig coords into sequence coords."""
    #  index 9 --> GC---CGC-T <-- index 18 (intervals are half closed)
    #              GCTAGCGCG- <-- the G is index 981
    self.assertEqual(self.pa1.alignment_to_sequence_coords("s1", 9, 19),
                     (107, 113))
    self.assertEqual(self.pa2.alignment_to_sequence_coords("s3", 9, 19),
                     (981, 990))
    #  index 12 --> --CGC <-- index 16 (intervals are half closed)
    self.assertEqual(self.pa1.alignment_to_sequence_coords("s1", 12, 17),
                     (109, 112))
    self.assertEqual(self.pa2.alignment_to_sequence_coords("s3", 12, 17),
                     (982, 987))
    #  index 9 --> GC-- <-- index 12 (intervals are half closed)
    self.assertEqual(self.pa1.alignment_to_sequence_coords("s1", 9, 13),
                     (107, 109))
    self.assertEqual(self.pa2.alignment_to_sequence_coords("s3", 9, 13),
                     (986, 990))
    #  index 24 --> --- <-- index 26 (intervals are half closed)
    self.assertEqual(self.pa1.alignment_to_sequence_coords("s2", 24, 27), None)
    self.assertEqual(self.pa2.alignment_to_sequence_coords("s3", 24, 27), None)
    # raise exception if invalid seq name
    self.assertRaises(InvalidSequenceError,
                      self.pa1.alignment_to_sequence_coords, "s3", 9, 19)
    # raise exception if start is greater than or equal to end
    self.assertRaises(InvalidAlignmentCoordinatesError,
                      self.pa1.alignment_to_sequence_coords, "s2", 15, 9)
    self.assertRaises(InvalidAlignmentCoordinatesError,
                      self.pa1.alignment_to_sequence_coords, "s2", 9, 9)
    # raise exception if coordinates are fuly outside alignemnt, ragrdless of
    # the trim option
    self.assertRaises(InvalidAlignmentCoordinatesError,
                      self.pa1.alignment_to_sequence_coords, "s2", 100, 150)
    self.assertRaises(InvalidAlignmentCoordinatesError,
                      self.pa1.alignment_to_sequence_coords, "s2", 100, 150,
                      trim=True)
    # raise exception when coordinates are prtially outside the alignment, but
    # not if the trim option is given.
    self.assertRaises(InvalidAlignmentCoordinatesError,
                      self.pa1.alignment_to_sequence_coords, "s1", 33, 40)
    r = self.pa1.alignment_to_sequence_coords("s1", 33, 40, trim=True)
    self.assertEqual(r, (126, 129))

  def test_liftover_coords(self):
    """Convert coords of one seq. in pairwise alig. into coords of other."""
    #  103 --> CGTAGC---CGC-T   <-- 112
    # 1014 --> CGTAGCTAGCGCG-   <-- 1016
    #  993 -->                  <-- 981
    self.assertEqual(self.pa1.liftover("s1", "s2", 103, 113),
                     [(1004, 1010), (1013, 1016)])
    self.assertEqual(self.pa2.liftover("s1", "s3", 103, 113),
                     [(988, 994), (982, 985)])
    self.assertEqual(self.pa1.liftover("s2", "s1", 1010, 1013), [])
    self.assertEqual(self.pa2.liftover("s3", "s1", 985, 988), [])
    # should fail if coordinates are outside the bounds of the sequence
    self.assertRaises(InvalidSequenceCoordinatesError, self.pa1.liftover,
                      "s1", "s2", 95, 113)
    # but should be okay if we allow them to be trimmed
    r = self.pa1.liftover("s1", "s2", 95, 113, trim=True)
    self.assertEqual(r, [(1001, 1010), (1013, 1016)])
    # unless it is entriely outside the sequence
    self.assertRaises(InvalidSequenceCoordinatesError, self.pa1.liftover,
                      "s1", "s2", 95, 99, trim=True)

  def test_failure_on_different_lengths(self):
    """Test failure when sequences passed to constructor are ragged."""
    s1 = Sequence("s1", "-TCGCGTAGC---CGC-TAGGATGCGAT-CTGA", 100, 127)
    s2 = Sequence("s2", "ATCGCGTAGCTAGCGCG-AGCTG---CGATGCT--", 1000, 1029)
    args = [s1, s2]
    self.assertRaises(MultipleAlignmentError, MultipleSequenceAlignment, args)

  def test_empty_seq_with_diff_length(self):
    """Test using set of sequences where an empty seq doesn't match in size."""
    s1 = Sequence("s1", "-TCGCGTAGC---CGC-TAGCTGATGCGAT-CTGA", 100, 129)
    s2 = Sequence("s2", "ATCGCGTAGCTAGCGCG-AGCTG---CGATGCT--", 1000, 1029)
    s3 = UnknownSequence("s3", 25049, 25049 + 1601, "+", 50103)
    msa = MultipleSequenceAlignment([s1, s2, s3])
    self.assertEqual(msa.get_column(1), {"s1": "-", "s2": "A"})

  def test_different_ungapped_legnths(self):
    """Test that we can build an MSA with seqs of diff ungapped length."""
    s1 = Sequence("s1", "------TCGCGTAGC", 100, 129)
    s2 = Sequence("s2", "---------CGCAGC", 1000, 1006)
    s3 = Sequence("s3", "ATCGCGT--------", 120, 127)
    msa = MultipleSequenceAlignment([s1, s2, s3])
    self.assertEqual(msa.get_column(7), {"s1": "T", "s2": "-", "s3": "T"})


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
