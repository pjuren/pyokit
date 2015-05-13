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
    if len(seq) != len(qual) :
      raise NGSReadError("failed to create FastqSequence object -- length " +
                         "of sequence data (" + str(len(seq)) + ")" +
                         "does not match length of quality string (" +
                         str(len(qual)) + ")")

    Sequence.__init__(self, name, seq, use_mut_str)

    self.seq_qual = qual

    # for quality scores
    self.LOWSET_SCORE = 64
    self.HIGHEST_SCORE = 104

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
    self.trimRight(len(self) - size)

  def trimRight(self, amount):
    """
      Trim this fastqSequence in-place by removing <amount> nucleotides from
      the 3' end (right end).

      :param amount: the number of nucleotides to trim from the right-side of
                     this sequence.
    """
    self.sequenceData = self.sequenceData[:-amount]
    self.sequenceQual = self.sequenceQual[:-amount]

  def trimLeft(self, amount):
    """
      Trim this fastqSequence in-place by removing <amount> nucleotides from
      the 5' end (left end).

      :param amount: the number of nucleotides to trim from the left-side of
                     this sequence.
    """
    self.sequenceData = self.sequenceData[amount:]
    self.sequenceQual = self.sequenceQual[amount:]

  def getRelativeQualityScore(self, i):
    val = self.sequenceQual[i]
    return (ord(val) - self.LOWSET_SCORE) / float(self.HIGHEST_SCORE -
                                                  self.LOWSET_SCORE)

  def reverseComplement(self, is_RNA=None):
    """
      Reverse complement this fastq sequence in-place.
    """
    Sequence.reverseComplement(self, is_RNA)
    self.sequenceQual = self.sequenceQual[::-1]

  def split(self, point=None):
    """
      Split this fastq sequence into two halves. The original sequence is left
      unaltered.

      :param point: the point (index) at which to split this sequence. If None
                    (the default), then we split in the middle.
      :return: two FastqSequence objects which correspond to the split of this
               sequence.
    """
    if point is None :
      point = len(self) / 2

    r1 = NGSRead(self.sequenceName + ".1", self.sequenceData[:point],
                 self.sequenceQual[:point])
    r2 = NGSRead(self.sequenceName + ".2", self.sequenceData[point:],
                 self.sequenceQual[point:])
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
    if self.sequenceName != other.sequenceName and not forceMerge :
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
    return "@" + self.sequenceName + "\n" + self.sequenceData +\
           "\n" + "+" + self.sequenceName + "\n" + self.sequenceQual


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

  def testeq(self):
    """
      test the equality operator for fastQ sequences.
    """
    r1 = NGSRead("ACTGCT", "s1", "BBBBBB")
    r2 = NGSRead("ACTGCT", "s1", "BBBBBB")
    r3 = NGSRead("ACTGCT", "s1", "BBBfBB")
    r4 = NGSRead("ACCGCT", "s1", "BBBBBB")
    r5 = NGSRead("TCTGCT", "s2", "fBBBBB")
    r6 = NGSRead("CCCCCC", "s3", "fBBBBB")
    r7 = NGSRead("CCCCCC", "s1", "fBBfBB")
    r7 = NGSRead("CCCCCC", "s6", "fBBfBB")

    self.assertTrue((r1 == r2) == True)   # same name, same seq, same qual
    self.assertTrue((r1 == r3) == False)  # same name, same seq, diff qual
    self.assertTrue((r1 == r4) == False)  # same name, diff seq, same qual
    self.assertTrue((r3 == r4) == False)  # same name, diff seq, diff qual
    self.assertTrue((r1 == r5) == False)  # diff name, diff seq, diff qual
    self.assertTrue((r5 == r6) == False)  # diff name, diff seq, same qual
    self.assertTrue((r6 == r7) == False)  # diff name, same seq, diff qual
    self.assertTrue((r3 == r4) == False)  # diff name, same seq, same qual

    self.assertTrue((r1 != r2) == False)  # same name, same seq, same qual
    self.assertTrue((r1 != r3) == True)   # same name, same seq, diff qual
    self.assertTrue((r1 != r4) == True)   # same name, diff seq, same qual
    self.assertTrue((r3 != r4) == True)   # same name, diff seq, diff qual
    self.assertTrue((r1 != r5) == True)   # diff name, diff seq, diff qual
    self.assertTrue((r5 != r6) == True)   # diff name, diff seq, same qual
    self.assertTrue((r6 != r7) == True)   # diff name, same seq, diff qual
    self.assertTrue((r3 != r4) == True)   # diff name, same seq, same qual

  def testClipadaptor(self):
    input_seq = NGSRead("ACTGCTAGCGATCGACT", "n1", "QQQQQQQQQQQQQQQQQ")
    adaptor = Sequence("adap", "AGCGATAGACT")
    expect = NGSRead("ACTGCTNNNNNNNNNNN", "n1", "QQQQQQQQQQQQQQQQQ")
    clip_adaptor(input_seq, adaptor)
    got = input_seq
    self.assertTrue(expect == got)


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
