#!/usr/bin/python

"""
  Date of Creation: 12th April 2015
  Description:   Defines class representing NGS reads

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
from pyokit.datastruct.sequence import Sequence


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class NGSReadError(Exception):
  """
  Class representing errors that occur when manipulating Fastq sequence objects
  """
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


###############################################################################
#                              NGS READ CLASS                                 #
###############################################################################

class NGSRead(Sequence):
  """
  Data structure for holding informaton about a next-generation sequencing
  read

  :param seq:         The nucleotide sequence data for the read, as a string.
                      Can be DNA or RNA. Note that there is no check to make
                      sure the sequence data is valid, that's the
                      responsibility of the caller.
  :param name:        A name describing the sequence. Can be any string.
  :param qual:        The quality string for this sequence -- must be the same
                      length as the nucleotide sequence.
  :param use_mut_str: Store the sequence data as a mutable string, rather than
                      a regular python string. This should make editing
                      operations much faster, but it comes at the expense of
                      less flexibility (e.g. the object can not be used as a
                      hash key because it is mutable.)

  :raise SequenceError: if the sequence data is not the same length as the
                        quality data.
  """

  def __init__(self, seq, name=None, qual=None, use_mut_str=False):
    """
      Constructor for FastqSequence class; see class level documentation for
      descriptions of parameters.
    """
    if len(seq) != len(qual):
      msg = ("failed to create FastqSequence object -- length of sequence " +
             "data (" + str(len(seq)) + ") does not match length of quality " +
             "string (" + str(len(qual)) + "). seq data: " + seq +
             "  qual data: " + qual)
      raise NGSReadError(msg)

    Sequence.__init__(self, name, seq, use_mut_str)

    self.seq_qual = qual

    # for quality scores
    # ILLUMINA 1.3+ Phred+64
    self.LOWSET_SCORE = 64
    self.HIGHEST_SCORE = 104
    # Illumina 1.8+ Phred+33
    self.LOWSET_SCORE_ILL_18_PHRD_33 = 33
    self.HIGHEST_SCORE_ILL_18_PHRD_33 = 74

  def __eq__(self, seq):
    """
      determine whether two fastqSequence objects are equal. To be equal, their
      sequence data (name, nuc. sequence) must match, as well as their quality
      data.

      :param seq: the other sequence to compare against.
      :return: True if this sequence is equal to seq, else False.
    """
    return (Sequence.__eq__(self, seq) and
            self.seq_qual == seq.seq_qual)

  def __ne__(self, read):
    """
      determine whether two fastqSequence objects are not equal. They are
      considered unequal if any of their sequence data (name, nuc. sequence)
      does not match, or if their quality data does not match.

      :param seq: the other sequence to compare against.
      :return: True if this sequence is not equal to seq, else False.
    """

    return (Sequence.__ne__(self, read) or
            self.seq_qual != read.seq_qual)

  def truncate(self, size):
    """
      truncate this fastqSequence in-place so it is only <size> nucleotides
      long

      :param size: the number of nucleotides to truncate to.
    """
    if size > len(self):
      raise NGSReadError("Trying to truncate NGS read to size " + str(size) +
                         ", but read is only " + str(len(self)) +
                         " nucleotides long")
    if size < 0:
      raise NGSReadError("Trying to truncate NGS read to size less than 0")

    self.trimRight(len(self) - size)

  def trimRight(self, amount):
    """
      Trim this fastqSequence in-place by removing <amount> nucleotides from
      the 3' end (right end).

      :param amount: the number of nucleotides to trim from the right-side of
                     this sequence.
    """
    if amount == 0:
      return
    self.sequenceData = self.sequenceData[:-amount]
    self.seq_qual = self.seq_qual[:-amount]

  def trimLeft(self, amount):
    """
      Trim this fastqSequence in-place by removing <amount> nucleotides from
      the 5' end (left end).

      :param amount: the number of nucleotides to trim from the left-side of
                     this sequence.
    """
    if amount == 0:
      return
    self.sequenceData = self.sequenceData[amount:]
    self.sequenceQual = self.sequenceQual[amount:]

  def getRelativeQualityScore(self, i, score_type="ILLUMINA_PHRED_PLUS_33"):
    """
    Get the realtive quality score (i.e. the phred quality score) for a given
    base.

    :raise: NGSReadError if the index is less than 0 or more than l-1 where l
            is the length of the read; NGSReadError if the encoded quality
            score at the given location is outside the expected range;
            NGSReadError if the read has no quality string; NGSReadError if
            sepficied encoding scheme is unknown.
    """
    # error out if no quality string, or index is out of range
    if self.seq_qual is None:
      raise NGSReadError("Error finding realtive quality for base call at " +
                         "location " + str(i) + " in read " + self.name +
                         "; read has no quality string")
    if i < 0 or i > self.seq_qual:
      raise NGSReadError("Error finding relative quality for base call at " +
                         "location " + str(i) + " in read " + self.name +
                         "; outside of range for read with quality score " +
                         "string of length " + str(len(self.seq_qual)))

    val = self.seq_qual[i]
    if score_type == "ILLUMINA_PHRED_PLUS_33":
      return ord(val) - self.LOWSET_SCORE_ILL_18_PHRD_33
    elif score_type == "ILLUMINA_PHRED_PLUS_64":
      return ord(val) - self.LOWSET_SCORE
    else:
      raise NGSReadError("Unknown score type: " + score_type)

  def reverse_complement(self, is_RNA=None):
    """
      Reverse complement this read in-place.
    """
    Sequence.reverseComplement(self, is_RNA)
    self.seq_qual = self.seq_qual[::-1]

  def split(self, point=None):
    """
    Split this read into two halves. Original sequence is left unaltered.

    The name of the resultant reads will have '.1' and '.2' appended to the
    name from the original read.

    :param point: the point (index, starting from 0) at which to split this
                  read -- everything before this index will be placed into the
                  first sequence, and everything at or after this index will be
                  placed in the second resultant sequence. If None
                  (the default), then we split in the middle; if the original
                  size is not a multiple of 2, the extra nucleotide is placed
                  into the second resultant sequence. Must be >= 0
                  and <= length of sequence.
    :return: two NGSRead objects which correspond to the split of this
             sequence.
    """
    if point is None:
      point = len(self) / 2
    if point < 0:
      raise NGSReadError("Cannot split read at index less than 0 " +
                         "(index provided: " + str(point) + ")")
    if point > len(self):
      raise NGSReadError("Cannot split read at index greater than read " +
                         "length (index provided: " + str(point) + ")")

    r1 = NGSRead(self.sequenceData[:point], self.name + ".1",
                 self.seq_qual[:point])
    r2 = NGSRead(self.sequenceData[point:], self.name + ".2",
                 self.seq_qual[point:])
    return r1, r2

  def merge(self, other, forceMerge=False):
    """
      Merge two reads by concatenating their sequence data and their
      quality data (<self> first, then <other>); <self> and <other> must have
      the same sequence name. A new merged FastqSequence object is returned;
      <Self> and <other> are left unaltered.

      :param other: the other sequence to merge with self.
      :param forceMerge: force the merge to occur, even if sequences names
                         don't match. In this case, <self> takes precedence.
      :return: A new FastqSequence that represents the merging of <self> and
               <other>
      :raise: FastqSequenceError if the sequences names do not match, and the
              forceMerge parameter is not set.
    """
    if self.sequenceName != other.sequenceName and not forceMerge:
      raise NGSReadError("cannot merge " + self.sequenceName + " with " +
                         other.sequenceName + " -- different " +
                         "sequence names")

    name = self.sequenceName
    seq = self.sequenceData + other.sequenceData
    qual = self.sequenceQual + other.sequenceQual

    return NGSReadError(name, seq, qual)

  def __str__(self):
    """
    :return: string representation of this NGS read
    """
    return self.to_fastq_str()

  def to_fastq_str(self):
    """
    :return: string representation of this NGS read in FastQ format
    """
    return "@" + self.name + "\n" + self.sequenceData +\
           "\n" + "+" + self.name + "\n" + self.seq_qual


###############################################################################
#                  HELPER FUNCTIONS FOR PROCESSING NGS READS                  #
###############################################################################

def clip_adaptor(read, adaptor):
  """
  Clip an adaptor sequence from this sequence. We assume it's in the 3'
  end. This is basically a convenience wrapper for clipThreePrime. It
  requires 8 out of 10 of the first bases in the adaptor sequence to match
  for clipping to occur.

  :param adaptor: sequence to look for. We only use the first 10 bases;
                  must be a full Sequence object, not just a string.
  """
  missmatches = 2
  adaptor = adaptor.truncate(10)
  read.clip_end(adaptor, len(adaptor) - missmatches)


def contains_adaptor(read, adaptor):
  """
  Check whether this sequence contains adaptor contamination. If it exists,
  we assume it's in the 3' end. This function requires 8 out of 10 of the
  first bases in the adaptor sequence to match for an occurrence to be
  reported.

  :param adaptor: sequence to look for. We only use first 10 bases; must be
                  a full Sequence object, not just string.
  :return: True if there is an occurence of <adaptor>, False otherwise
  """
  origSeq = read.sequenceData
  clip_adaptor(read, adaptor)
  res = False
  if read.sequenceData != origSeq:
    res = True
  read.sequenceData = origSeq
  return res


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class NGSReadUnitTests(unittest.TestCase):
  """
    Unit tests for NGSRead module
  """
  def setUp(self):
    self.r1 = NGSRead("ACTGCT", "s1", "BBBBBB")
    self.r2 = NGSRead("ACTGCT", "s1", "BBBBBB")
    self.r3 = NGSRead("ACTGCT", "s1", "BBBfBB")
    self.r4 = NGSRead("ACCGCT", "s1", "BBBBBB")
    self.r5 = NGSRead("TCTGCT", "s2", "fBBBBB")
    self.r6 = NGSRead("CCCCCC", "s3", "fBBBBB")
    self.r7 = NGSRead("CCCCCC", "s1", "fBBfBB")

  def test_length_mismatch(self):
    """Test trying to make a read with mismatched seq and qual data."""
    self.assertRaises(NGSReadError, NGSRead, "ACTGCT", "s1", "BBBBB")
    self.assertRaises(NGSReadError, NGSRead, "ACTGC", "s1", "BBBBBB")
    self.assertRaises(NGSReadError, NGSRead, "ACTGCT", "s1", "")
    self.assertRaises(NGSReadError, NGSRead, "", "s1", "BBBBB")

  def test_truncate(self):
    """Test truncating NGSReads."""
    self.r1.truncate(1)
    self.r2.truncate(6)
    self.r3.truncate(0)

    self.assertEqual(self.r1, NGSRead("A", "s1", "B"))
    self.assertEqual(self.r2, NGSRead("ACTGCT", "s1", "BBBBBB"))
    self.assertEqual(self.r3, NGSRead("", "s1", ""))
    self.assertRaises(NGSReadError, self.r4.truncate, -1)
    self.assertRaises(NGSReadError, self.r5.truncate, 7)

  def test_rev_comp(self):
    self.r1.reverse_complement()
    self.r3.reverse_complement()
    self.r7.reverse_complement()

    self.assertEqual(self.r1, NGSRead("AGCAGT", "s1", "BBBBBB"))
    self.assertEqual(self.r3, NGSRead("AGCAGT", "s1", "BBfBBB"))
    self.assertEqual(self.r7, NGSRead("GGGGGG", "s1", "BBfBBf"))

  def test_split(self):
    self.assertEquals(self.r1.split(), (NGSRead("ACT", "s1.1", "BBB"),
                                        NGSRead("GCT", "s1.2", "BBB")))
    self.assertEquals(self.r1.split(3), (NGSRead("ACT", "s1.1", "BBB"),
                                         NGSRead("GCT", "s1.2", "BBB")))
    self.assertEquals(self.r1.split(0), (NGSRead("", "s1.1", ""),
                                         NGSRead("ACTGCT", "s1.2", "BBBBBB")))
    self.assertEquals(self.r1.split(6), (NGSRead("ACTGCT", "s1.1", "BBBBBB"),
                                         NGSRead("", "s1.2", "")))
    self.assertRaises(NGSReadError, self.r1.split, -1)
    self.assertRaises(NGSReadError, self.r1.split, 7)

  def testeq(self):
    """
      test the equality operator for NGSRead sequences.
    """
    self.assertTrue((self.r1 == self.r2) is True)   # same name, seq, qual
    self.assertTrue((self.r1 == self.r3) is False)  # same name, seq, diff qual
    self.assertTrue((self.r1 == self.r4) is False)  # same name, qual, diff seq
    self.assertTrue((self.r3 == self.r4) is False)  # same name, diff seq, qual
    self.assertTrue((self.r1 == self.r5) is False)  # diff name, seq, qual
    self.assertTrue((self.r5 == self.r6) is False)  # diff name, seq, same qual
    self.assertTrue((self.r6 == self.r7) is False)  # diff name, qual, same seq
    self.assertTrue((self.r3 == self.r4) is False)  # diff name, same seq, qual

    self.assertTrue((self.r1 != self.r2) is False)  # same name, seq, same qual
    self.assertTrue((self.r1 != self.r3) is True)   # same name, seq, diff qual
    self.assertTrue((self.r1 != self.r4) is True)   # same name, qual, diff seq
    self.assertTrue((self.r3 != self.r4) is True)   # same name, diff seq, qual
    self.assertTrue((self.r1 != self.r5) is True)   # diff name, diff seq, qual
    self.assertTrue((self.r5 != self.r6) is True)   # diff name, seq, same qual
    self.assertTrue((self.r6 != self.r7) is True)   # diff name, qual, same seq
    self.assertTrue((self.r3 != self.r4) is True)   # diff name, same seq, qual

  def testClipadaptor(self):
    input_seq = NGSRead("ACTGCTAGCGATCGACT", "n1", "QQQQQQQQQQQQQQQQQ")
    adaptor = Sequence("adap", "AGCGATAGACT")
    expect = NGSRead("ACTGCTNNNNNNNNNNN", "n1", "QQQQQQQQQQQQQQQQQ")
    clip_adaptor(input_seq, adaptor)
    got = input_seq
    self.assertTrue(expect == got)

  def testGetRelativeQualityScore(self):
    input_seq = NGSRead("TTTGTAGTTTTGATTTTTGAGGTTTTAGTTATTTATTAGTAAATTAAGGT" +
                        "TTAAAAATATTAAATGGAAAATTTTAGAAATAAGTAATGTATAAGTTTTAA",
                        "SN346:615:HMTFNADXX:1:1101:1413:2147 1:N:0:CGATGT",
                        "BBBFFFFFFFFFFFIIIIFFFIIFFIFFFFIIIIIIIIIIIFIIIIFFII" +
                        "IIFIIIFFFFIFFFFFFFFFFFFFIIFIFIFFFFFFFFFFBBFBBFBFFBB")
    self.assertEquals(input_seq.getRelativeQualityScore(0), 33)
    self.assertEquals(input_seq.getRelativeQualityScore(3), 37)


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
