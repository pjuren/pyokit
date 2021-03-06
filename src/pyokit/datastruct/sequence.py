#!/usr/bin/python

"""
  Date of Creation: 13th September 2010
  Description:   Defines superclass for fastq and fasta reads

  Copyright (C) 2010-2014
  Philip J. Uren

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

import unittest
from pyokit.util.progressIndicator import ProgressIndicator
from pyokit.util.mutableString import MutableString

###############################################################################
#                             MODULE CONSTANTS                                #
###############################################################################

DNA_COMPLEMENTS = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
                   "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}
RNA_COMPLEMENTS = {"A": "U", "U": "A", "C": "G", "G": "C", "N": "N",
                   "a": "u", "u": "a", "c": "g", "g": "c", "n": "n"}
GAP_CHAR = "-"
UNKNOWN_SEQ_NAME = "UNKNOWN_SEQUENCE"
RNA_NUCS = "ACGUNacgun"
DNA_NUCS = "ACGTNacgtn"


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class SequenceError(Exception):
  """
  Class representing errors that occur when manipulating general sequence
  objects
  """
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


class InvalidSequenceCoordinatesError(SequenceError):
  """
  Thrown when invalid values are given to index into a sequence.
  """
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


###############################################################################
#                               SEQUENCE CLASS                                #
###############################################################################

class Sequence(object):
  """
    This is the base class for all sequences in Pyokit. Objects from this class
    will have only a sequence name and actual nucleotide sequence data.

    :param seqName:          A name describing the sequence. Can be any string.
    :param seqData:          The nucleotide sequence data. Can be DNA or RNA.
                             Note that there is no check to make sure the
                             sequence data is valid, that's the responsibility
                             of the caller.
    :param start_coord:      TODO
    :param end_coord:        TODO
    :param strand:           By default, this is +, but can also be set to -
                             to indicate that this sequence is a reverse
                             complement.
    :param remaining:        the amount of sequence that comes after this; 0 if
                             this is the whole sequence. Alterntively, you
                             might think of this as the negative strand
                             coordinates of the end of this sequence.
    :param meta_data:        dictionary containing meta-data key-value pairs
    :param useMutableString: Store the sequence data as a mutable string,
                             rather than a regular python string. This should
                             make editing operations must faster, but it comes
                             at the expense of less flexibility (e.g. the
                             object can not be used as a hash key because it
                             is mutable.)
  """

  def __init__(self, seqName, seqData, start_coord=None, end_coord=None,
               strand="+", remaining=0, meta_data=None,
               useMutableString=False):
    """
      Constructor for Sequence objects. See class level documentation for
      parameter descriptions.
    """
    if strand != "+" and strand != "-":
      raise ValueError("Sequence strand must be either + or -, found " +
                       strand + " instead")
    self.name = seqName
    if useMutableString:
      self.sequenceData = MutableString(seqData)
    else:
      self.sequenceData = seqData
    self.mutableString = useMutableString
    self.remaining = remaining
    self._start_coord = start_coord
    self._end_coord = end_coord
    self.strand = strand
    self._ungapped_len = None   # we compute this just-in-time..
    self._effective_len = None  # .. and this
    self.meta_data = meta_data if meta_data is not None else {}

  def copy(self):
    """
    Copy constructor for Sequence objects.
    """
    return Sequence(self.name, self.sequenceData, self.start, self.end,
                    self.strand, self.remaining, self.meta_data,
                    self.mutableString)

  @property
  def start(self):
    """
    :return: The coordinate of the first nucleotide in this sequence; by
    convention, we call this coordinate 1 if no other value was provided.
    """
    if self._start_coord is None:
      return 1
    return self._start_coord

  @property
  def end(self):
    """
    :return: The coordinate of the end of this sequence; as with all other
    indexing of sequences in pyokit, sequences are not inclusive of their last
    index. Computed just-in-time from the ungapped sequence length if it
    wasn't provided at construction time.
    """
    if self._end_coord is None:
      return self.ungapped_len + 1
    return self._end_coord

  @property
  def ungapped_len(self):
    if self._ungapped_len is None:
      self._ungapped_len = 0
      for nuc in self.sequenceData:
        if nuc != GAP_CHAR:
          self._ungapped_len += 1
    # take this oportunity to check that coords match ungapped sequence len
    e_ok = self._end_coord is not None
    if e_ok and self._ungapped_len != self.end - self.start:
      msg = ("ungapped length (" + str(self._ungapped_len) + ") of sequence " +
             "doesn't match start (" + str(self.start) + ") and end (" +
             str(self.end) + ") coordinates")
      raise SequenceError(msg)
    return self._ungapped_len

  @property
  def effective_len(self):
    """
    Get the length of the sequence if N's are disregarded.
    """
    if self._effective_len is None:
      self._effective_len = len([nuc for nuc in self.sequenceData
                                 if nuc != "N" and nuc != "n"])
    return self._effective_len

  def __len__(self):
    """
    Get the length of the sequence, defined as the length of its sequence
    data
    """
    return len(self.sequenceData)

  def __getitem__(self, i):
    return self.sequenceData[i]

  def subsequence(self, start, end):
    """
    Extract a subsequence from this sequence object using absolute coordinates
    that exist in the same coordinate space as the sequence itself. For
    example:

      46  --> A--CTGC-TAGC-GATCGACT <--  62      subsequence(47,52) == CTGC-T

    :param start:      the index marking the start (inclusive) of the
                       subsequence. This is a one-based index, and is in the
                       same coordinate space as this sequence object.
    :param end:        the index marking the end (exclusive) of the
                       subsequence. This is a one-based index, and is in the
                       same coordinate space as this sequence object.

    :return: a new sequence object that represents the subsequence of this from
             position start (indexed from 1, inclusive) to end (indexed from
             1, exclusive). Ungapped length of this will always be equal to
             end - start
    :rasie SequenceError: if the coordinates given fall outside of the start
                          and end indices of this sequence object.
    """
    if (start < self.start or start >= self.end or end <= self.start or
       end > self.end or start >= end):
      raise SequenceError("invalid subsequence coordinates (" + str(start) +
                          " ," + str(end) + ") for sequence " + str(self))

    rel_start_string_coord = self.start - start - 1
    rel_end_string_coord = self.end - end - 1
    seq = self.sequenceData[rel_start_string_coord:rel_end_string_coord]
    return Sequence(self.name, seq, start, end, self.strand,
                    self.remaining + self.end - end, self.mutableString)

  def relative_subsequence(self, start, end):
    """
    Extract a subsequence from this sequence using coordinates that are
    relative to the start (relative position 1) and end coordinates of the
    sequence. For example:

      46  --> A--CTGC-TAGC-GATCGACT <--  62      subsequence(2,7) == CTGC-T

    :param start:      the index marking the start (inclusive) of the
                       subsequence. This is a one-based index, and is in the
                       coordinate space of the sequence (i.e. from 1 to N,
                       where N is the number of non-gap nucleotides in the
                       sequence)
    :param end:        the index marking the end (exclusive) of the
                       subsequence. This is a one-based index, and is in the
                       coordinate space of the sequence (i.e. from 1 to N,
                       where N is the number of non-gap nucleotides in the
                       sequence)

    :return: a new sequence object that represents the subsequence
    :rasie SequenceError: if the start coordinate is less than 1 or the end
                          coordinate is greater than the ungapped length of
                          this sequence.
    """
    if start < 1:
      raise SequenceError("invalid start coordinate for subsequence: " +
                          str(start))
    if end > self.ungapped_len:
      raise SequenceError("invalid end coordinate for subsequence: " +
                          str(end) + " greater than number of non-gap " +
                          "nucleotides in sequence (" + self.ungapped_len +
                          ")")
    raise SequenceError("method not implemented")

  def gapped_relative_subsequence(self, start, end):
    """
    Extract a subsequence from this sequence using coordinates that are
    relative to the start of the sequence (relative position 1) and the number
    of nuceltodies in the sequence, including gaps. For example:

      46  --> A--CTGC-TAGC-GATCGACT <--  62      subsequence(2,7) == --CTG
    """
    if start < 1:
      msg = "invalid start coordinate for subsequence: " + str(start)
      raise InvalidSequenceCoordinatesError(msg)
    if end > len(self.sequenceData) + 1:
      msg = "invalid end coordinate for subsequence: " + str(end) +\
            " greater than length of sequence (" + str(len(self)) + ")"
      raise InvalidSequenceCoordinatesError(msg)

    non_gaps_before = start - 1 - self.sequenceData[:start - 1].count(GAP_CHAR)
    seq = self.sequenceData[start - 1:end - 1]
    non_gaps_in = len(seq) - seq.count(GAP_CHAR)
    non_gaps_after = self.ungapped_len - non_gaps_before - non_gaps_in
    new_start = self.start + non_gaps_before
    new_end = new_start + non_gaps_in
    new_remaining = self.remaining + non_gaps_after

    return Sequence(self.name, seq, new_start, new_end, self.strand,
                    new_remaining, self.meta_data, self.mutableString)

  def is_positive_strand(self):
    """
    """
    return self.strand == "+"

  def percentNuc(self, nuc):
    """
      return the percentage of the sequence which is equal to the passed nuc.

      :param nuc: the nucleotide to compute percentage composition for. There
                  is no check to make sure this is a valid nucleotide.
      :return: the percentage of the sequence that is <nuc>
    """
    count = reduce(lambda x, y: x + 1 if y == nuc else x, self.sequenceData, 0)
    return count / float(len(self.sequenceData))

  def similarity(self, self_start, self_end, other_start, other_end, other):
    """
      Compute the number of matching bases in the subsequences self[start, end]
      and other[o_start, o_end]. Note that the subsequences must be the same
      length.

      :param self_start:  start index for sub-sequence in self
      :param self_end:    end index for sub-sequence in self
      :param other_start: start index for subsequence in other sequence
      :param other_end:   end index for subsequence in other sequence
      :param other:       other sequence to compare to this.
    """
    assert(self_end - self_start == other_end - other_start)
    count = 0
    for i in range(0, self_end - self_start + 1):
      if (self.sequenceData[self_start + i] ==
         other.sequenceData[other_start + i]):
        count += 1
    return count

  def reverseComplement(self, isRNA=None):
    """
      Reverse complement this sequence in-place.

      :param isRNA: if True, treat this sequence as RNA. If False, treat it as
                    DNA. If None (default), inspect the sequence and make a
                    guess as to whether it is RNA or DNA.
    """
    isRNA_l = self.isRNA() if isRNA is None else isRNA

    tmp = ""
    for n in self.sequenceData:
      if isRNA_l:
        tmp += RNA_COMPLEMENTS[n]
      else:
        tmp += DNA_COMPLEMENTS[n]
    self.sequenceData = tmp[::-1]

  def __eq__(self, seq):
    """
      Check wheter this sequence is equal to another sequence. Sequences are
      equal if they have the same name, nucleotide sequence, coordinates,
      strand and meta data.

      :param seq: the other sequence to compare against.
      :return: true if this sequence is equal to passed parameter, else false.
    """
    if seq is None:
      return False
    return (self.sequenceData == seq.sequenceData and
            self.name == seq.name and
            self.meta_data == seq.meta_data and
            self.start == seq.start and
            self.end == seq.end and
            self.remaining == seq.remaining)

  def __ne__(self, read):
    """
      Check wheter this sequence is not equal to another sequence. Sequences
      are equal if they have the same name and nucleotide sequence.

      :param seq: the other sequence to compare against.
      :return: true if this sequence is not equal to passed param., else false.
    """
    if read is None:
      return True
    return (self.sequenceData != read.sequenceData or
            self.name != read.name)

  def nsLeft(self, amount):
    """
      Replace leftmost <amount> bases by Ns.
    """
    self.sequenceData = (amount * "N") + self.sequenceData[amount:]

  def nsRight(self, amount):
    """
      Replace rightmost <amount> bases by Ns
    """
    self.sequenceData = self.sequenceData[:-amount] + (amount * "N")

  def maskRegion(self, region):
    """
      Replace nucleotides in this sequence in the regions given by Ns

      :param region: any object with .start and .end attributes. Co-ords are
                     zero based and inclusive of both end points. Any other
                     attributes (e.g. chrom.) are ignored.
      :raise SequenceError: if region specifies nucleotides not present in
                            this sequence
    """
    if region.start < 0 or region.end < 0 or \
       region.start > len(self) or region.end > len(self):
      raise SequenceError("cannot mask region " + str(region.start) + " to " +
                          str(region.end) + " in " + self.name + ". " +
                          "Region specifies nucleotides not present in " +
                          "this read. Valid range would have been 0 -- " +
                          str(len(self)))

    if self.mutableString:
      for i in range(region.start, region.end + 1):
        self.sequenceData[i] = 'N'
    else:
      self.sequenceData = "".join([self.sequenceData[:region.start],
                                  ("N" * (region.end - region.start + 1)),
                                  self.sequenceData[region.end + 1:]])

  def maskRegions(self, regions, verbose=False):
    """
      Mask the given regions in this sequence with Ns.

      :param region: iterable of regions to mask. Each region can be any object
                     with .start and .end attributes. Co-ords are zero based
                     and inclusive of both end points. Any other attributes
                     (e.g. chrom.) are ignored.
      :param verbose: print status messages to stderr if True
    """
    if verbose:
      pind = ProgressIndicator(totalToDo=len(regions),
                               messagePrefix="completed",
                               messageSuffix="of masking regions in " +
                                             self.name)
    for region in regions:
      self.maskRegion(region)
      if verbose:
        pind.done += 1
        pind.showProgress()

  def isDNA(self):
    """
      Make a guess as to whether this sequence is a DNA sequence or not by
      looking at the symbols it contains.

      :return: True if contains only DNA nucleotides, False otherwise
    """
    for nuc in self.sequenceData:
      if nuc not in DNA_NUCS:
        return False
    return True

  def isRNA(self):
    """
      Make a guess as to whether this sequence is an RNA sequence or not by
      looking at the symbols it contains.

      :return: True if contains only RNA nucleotides, False otherwise
    """
    for nuc in self.sequenceData:
      if nuc not in RNA_NUCS:
        return False
    return True

  def toRNA(self):
    """
      Convert this sequence in-place to an RNA sequence by changing any Ts
      to Us
    """
    self.sequenceData = self.sequenceData.replace("T", "U")

  def toDNA(self):
    """
      Convert this sequence in-place to a DNA sequence by changing any Us to Ts
    """
    self.sequenceData = self.sequenceData.replace("U", "T")

  def split(self, point=None):
    """
      Split this sequence into two halves and return them. The original
      sequence remains unmodified.

      :param point: defines the split point, if None then the centre is used
      :return: two Sequence objects -- one for each side
    """
    if point is None:
      point = len(self) / 2

    r1 = Sequence(self.name + ".1", self.sequenceData[:point])
    r2 = Sequence(self.name + ".2", self.sequenceData[point:])

    return r1, r2

  def truncate(self, newLength):
    """
      Truncate this sequence in-place so it's only <newLength> nucleotides
      long.

      :param newLength: the length to truncate this sequence to.
    """
    return Sequence(self.name, self.sequenceData[:newLength])

  def clip_end(self, seq, mm_score):
    """
      Clip a sequence from the end of this sequence -- we assume the sequence
      to be clipped will always begin somewhere in this sequence, but may not
      be fully contained. If found, replaced with Ns.

      :param seq: sequence to be clipped
      :param mm_score: the number of matching bases needed to consider a hit,
                       mm_score = len(seq) would be 100% match
    """
    lim = mm_score - 1
    other_end = len(seq) - 1
    other_start = 0

    for i in range(len(self.sequenceData) - 1, lim - 1, -1):
      self_end = i
      self_start = i - (len(seq) - 1)

      if self_start < 0:
        self_start = 0
        other_start = other_end - self_end

      score = self.similarity(self_start, self_end, other_start,
                              other_end, seq)
      if (score >= mm_score):
        self.nsRight(len(seq) + (len(self) - i) - 1)
        break

  def isPolyA(self):
    """
      Determine whether this sequence is polyA. To be a polyA sequence, it
      must have > 90% Adenine.

      :return: True if the sequence is PolyA by the above definition.
    """
    return self.percentNuc("A") >= 0.9

  def isPolyT(self):
    """
      Determine whether this sequence is polyT. To be a polyT sequence, it
      must have > 90% Thymine.

      :return: True if the sequence is PolyT by the above definition.
    """
    return self.percentNuc("T") >= 0.9

  def isLowQuality(self):
    """
      Determine whether this is a low quality sequence. To be considered a low
      quality sequence, it must have > 10% Ns.

      :return: True if this sequence meets the above definition of low-quality.
    """
    return self.percentNuc("N") >= 0.1

  def maskMatch(self, mask):
    """
      Determine whether this sequence matches the given mask.

      :param mask: string to match against. Ns in the mask are considered to
                   match anything in the sequence -- all other chars must
                   match exactly.
      :return: True if the mask matches at all places, otherwise false
    """
    if len(mask) > len(self.sequenceData):
      return False
    lim = len(mask)
    for i in range(0, lim):
      if mask[i] == "N" or mask[i] == "n":
        continue
      if mask[i] != self.sequenceData[i]:
        return False
    return True

  #############################################################################
  #                     SEQUENCE CLASS -- STRING FORMATTING                   #
  #############################################################################

  def __str__(self):
    """
    :return: string representation of this sequence object
    """
    return self.to_fasta_str()

  def meta_data_to_string(self):
    """
    """
    res = ""
    first = True
    for k in self.meta_data:
      if first:
        first = False
      else:
        res += " "
      res += (str(k) + "=" + str(self.meta_data[k]))
    return res

  def to_fasta_str(self, line_width=50, include_coords=True):
    """
    :return: string representation of this sequence object in fasta format
    """
    res = ">" + self.name
    if include_coords:
      res += ":" + str(self.start) + "-" + str(self.end)
      res += " (" + str(self.remaining) + ")"
    m_str = self.meta_data_to_string()
    if m_str:
      res += (" " + m_str)
    res += "\n"
    for i in range(0, len(self.sequenceData), line_width):
      res += self.sequenceData[i:i + line_width]
      if i + line_width < len(self.sequenceData):
        res += "\n"
    return res


###############################################################################
#                              UNKNOWN SEQUENCE                               #
###############################################################################
class UnknownSequence(Sequence):
  """
  Represents a sequence where meta data about the sequence is known, but the
  actual sequence data is unknown. Behaves much like a regular Sequence, but
  any attempts to access the actual sequence data results in exceptions.
  """
  def __init__(self, seqName, start_coord=None, end_coord=None,
               strand="+", remaining=0, meta_data=None,
               useMutableString=False):
    Sequence.__init__(self, seqName, "", start_coord, end_coord, strand,
                      remaining, meta_data, useMutableString)

  @property
  def ungapped_len(self):
    return len(self)

  @property
  def effective_len(self):
    """
    Get the length of the sequence if N's are disregarded.
    """
    return len(self)

  def __len__(self):
    """
    Get the length of the sequence, defined as the length of its sequence
    data
    """
    return self.end - self.start

  def __getitem__(self, i):
    raise SequenceError("No sequence data")

  def subsequence(self, start, end):
    raise SequenceError("No sequence data")

  def relative_subsequence(self, start, end):
    raise SequenceError("No sequence data")

  def gapped_relative_subsequence(self, start, end):
    raise SequenceError("No sequence data")


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class SequenceUnitTests(unittest.TestCase):
  """
    Unit tests for sequence class
  """

  def test_gapped_relative_subsequence(self):
    """
      46 -->  A--CTGC-TAGC-GATCGACT  <-- 62   **  subsequence(2,7) == --CTG

      def __init__(self, seqName, seqData, start_coord=None, end_coord=None,
                   strand="+", remaining=0, useMutableString=False):
    """
    s1 = Sequence("s1", "A--CTGC-TAGC-GATCGACT", 46, 63, "+", 10)
    r1 = s1.gapped_relative_subsequence(2, 7)
    self.assertEqual(r1, Sequence("s1", "--CTG", 47, 50, "+", 23))
    self.assertEqual(r1.ungapped_len, 3)

  def testNsLeft(self):
    input_seq = Sequence("name", "ACTGCTAGCGATCGACT")
    expect = Sequence("name", "NNNNNTAGCGATCGACT")
    input_seq.nsLeft(5)
    got = input_seq
    self.assertTrue(expect == got)

  def testNsRight(self):
    input_seq = Sequence("name", "ACTGCTAGCGATCGACT")
    expect = Sequence("name", "ACTGCTAGCGATNNNNN")
    input_seq.nsRight(5)
    got = input_seq
    self.assertTrue(expect == got)

  def testLengths(self):
    input_seq = Sequence("name", "ACTNCTANCGATNNACT")
    self.assertTrue(len(input_seq) == 17)
    self.assertTrue(input_seq.effective_len == 13)

  def testMaskRegion(self):
    class TestRegion:
      def __init__(self, s, e):
        self.start = s
        self.end = e

    input_seq = Sequence("name", "ACTNCTANCGATNNACT")
    expect = Sequence("name", "ANNNNTANCGATNNACT")
    out = input_seq.copy()
    out.maskRegion(TestRegion(1, 4))
    self.assertTrue(expect == out)

  def testCopyConstructor(self):
    one = Sequence("name", "ACTNCTANCGATNNACT")
    two = one.copy()
    self.assertTrue(one == two)
    one.name = "old"
    self.assertTrue(one != two)

  def testReverseComplement(self):
    input_seq = Sequence("name", "ACTGCTAGCATGCGNN")
    expect = Sequence("name", "NNCGCATGCTAGCAGT")
    input_seq.reverseComplement()
    self.assertTrue(input_seq == expect)

  def testFormattedString(self):
    """
      test that string formatting works correctly for fasta sequences
    """
    r = Sequence("name", "ATCGATCGATCGATCTCGA")
    expect = ">name:1-20 (0)\n" +\
             "ATCGA\n" +\
             "TCGAT\n" +\
             "CGATC\n" +\
             "TCGA"
    got = r.to_fasta_str(line_width=5)
    self.assertEqual(got, expect)

    # make sure this also works with a mutable underlying sequence
    r = Sequence("name", "ATCGATCGATCGATCTCGA", useMutableString=True)
    expect = ">name:1-20 (0)\n" +\
             "ATCGA\n" +\
             "TCGAT\n" +\
             "CGATC\n" +\
             "TCGA"
    got = r.to_fasta_str(line_width=5)
    self.assertEqual(got, expect)


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
