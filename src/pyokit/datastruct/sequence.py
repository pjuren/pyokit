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

import unittest, sys, os
from pyokit.util.progressIndicator import ProgressIndicator

class MutableString :
  """
    Strings in python are immutable. That brings a number of advantages, but
    one problem is that they are expensive to edit. This class implements a
    string as a list of char, which is cheaper to edit, but cannot be used for
    things like dictionary keys due to it's mutability.
  """
  def __init__(self, strng):
    self.list = []
    for ch in strng :
      self.list.append(ch)

  def __getitem__(self, indx):
    if isinstance(indx, slice):
      return "".join(self.list.__getitem__(indx))
    return self.list[indx]

  def __setitem__(self, k, v):
    self.list[k] = v

  def __len__(self):
    return len(self.list)

  def __eq__(self, other):
    if type(other).__name__ == "str" :
      other = [x for x in other]
    else :
      try :
        other = other.list
      except : return false
    return self.list == other

  def __ne__(self, other):
    if type(other).__name__ == "str" :
      other = [x for x in other]
    else :
      try :
        other = other.list
      except : return false
    return self.list != other

  def __str__(self):
    return "".join(self.list)


class Sequence:
  """
    This is the base class for all sequences in Pyokit. Objects from this class
    will have only a sequence name and actual nucleotide sequence data.

    :param seqName:          A name describing the sequence. Can be any string.
    :param seqData:          The nucleotide sequence data. Can be DNA or RNA.
                             Note that there is no check to make sure the
                             sequence data is valid, that's the responsibility
                             of the caller.
    :param useMutableString: Store the sequence data as a mutable string, rather
                             than a regular python string. This should make
                             editing operations must faster, but it comes at the
                             expense of less flexibility (e.g. the object can
                             not be used as a hash key because it is mutable.)
  """


  DNA_COMPLEMENTS = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N",
                     "a":"t", "t":"a", "c":"g", "g":"c", "n":"n"}
  RNA_COMPLEMENTS = {"A":"U", "U":"A", "C":"G", "G":"C", "N":"N",
                     "a":"u", "u":"a", "c":"g", "g":"c", "n":"n"}

  def __init__(self, seqName, seqData, useMutableString = False):
    """
      Constructor for Sequence objects. See class level documentation for
      parameter descriptions.
    """
    self.sequenceName = seqName
    if useMutableString :
      self.sequenceData = MutableString(seqData)
    else :
      self.sequenceData = seqData
    self.mutableString = useMutableString

  def copy(self):
    """
      Copy constructor for Sequence objects.
    """
    return Sequence(self.sequenceName, self.sequenceData, self.mutableString)

  def percentNuc(self, nuc):
    """
      return the percentage of the sequence which is equal to the passed nuc.

      :param nuc: the nucleotide to compute percentage composition for. There is
                  no check to make sure this is a valid nucleotide.
      :return: the percentage of the sequence that is <nuc>
    """
    count = reduce(lambda x,y: x+1 if y==nuc else x, self.sequenceData, 0)
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
    for i in range(0, self_end - self_start + 1) :
      if self.sequenceData[self_start + i] == other.sequenceData[other_start + i] :
        count += 1
    return count

  def reverseComplement(self, isRNA=None):
    """
      Reverse complement this sequence in-place.

      :param isRNA: if True, treat this sequence as RNA. If False, treat it as
                    DNA. If None (default), inspect the sequence and make a
                    guess as to whether it is RNA or DNA.
    """
    isRNA_l = self.isRNA() if isRNA == None else isRNA

    tmp = ""
    for n in self.sequenceData :
      if isRNA_l : tmp += Sequence.RNA_COMPLEMENTS[n]
      else : tmp += Sequence.DNA_COMPLEMENTS[n]
    self.sequenceData = tmp[::-1]

  def __len__(self):
    """
      Get the length of the sequence, defined as the length of its sequence data
    """
    return len(self.sequenceData)

  def effectiveLength(self):
    """
      Get the length of the sequence if N's are disregarded.
    """
    return len([nuc for nuc in self.sequenceData
                    if nuc != "N" and nuc != "n"])

  def __eq__(self, seq):
    """
      Check wheter this sequence is equal to another sequence. Sequences are
      equal if they have the same name and nucleotide sequence.

      :param seq: the other sequence to compare against.
      :return: true if this sequence is equal to passed parameter, else false.
    """
    if seq == None : return False
    return  self.sequenceData == seq.sequenceData and\
            self.sequenceName == seq.sequenceName

  def __ne__(self, read):
    """
      Check wheter this sequence is not equal to another sequence. Sequences are
      equal if they have the same name and nucleotide sequence.

      :param seq: the other sequence to compare against.
      :return: true if this sequence is not equal to passed param., else false.
    """
    if read == None : return True
    return  self.sequenceData != read.sequenceData or\
            self.sequenceName != read.sequenceName

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
       region.start > len(self) or region.end > len(self) :
      raise SequenceError("cannot mask region " + str(region.start) + " to " +\
                          str(region.end) + " in " + self.sequenceName + ". " +\
                          "Region specifies nucleotides not present in " +\
                          "this read. Valid range would have been 0 -- " +\
                          str(len(self)))

    if self.mutableString :
      for i in range(region.start, region.end + 1) : self.sequenceData[i] = 'N'
    else :
      self.sequenceData = "".join([self.sequenceData[:region.start],
                          ("N" * (region.end - region.start + 1)),
                          self.sequenceData[region.end+1:]])

  def maskRegions(self, regions, verbose = False):
    """
      Mask the given regions in this sequence with Ns.

      :param region: iterable of regions to mask. Each region can be any object
                     with .start and .end attributes. Co-ords are zero based and
                     inclusive of both end points. Any other attributes (e.g.
                     chrom.) are ignored.
      :param verbose: print status messages to stderr if True
    """
    if verbose:
      pind = ProgressIndicator(totalToDo = len(regions),
                               messagePrefix = "completed",
                               messageSuffix = "of masking regions in " +\
                                                self.sequenceName)
    for region in regions :
      self.maskRegion(region)
      if verbose :
        pind.done += 1
        pind.showProgress()

  def isDNA(self):
    """
      Make a guess as to whether this sequence is a DNA sequence or not by
      looking at the symbols it contains.

      :return: True if contains only DNA nucleotides, False otherwise
    """
    for nuc in self.sequenceData :
      if not nuc in "ACGTacgtn" : return False
    return True

  def isRNA(self):
    """
      Make a guess as to whether this sequence is an RNA sequence or not by
      looking at the symbols it contains.

      :return: True if contains only RNA nucleotides, False otherwise
    """
    for nuc in self.sequenceData :
      if not nuc in "ACGUacgun" : return False
    return True

  def toRNA(self):
    """
      Convert this sequence in-place to an RNA sequence by changing any Ts to Us
    """
    self.sequenceData = self.sequenceData.replace("T","U")

  def toDNA(self):
    """
      Convert this sequence in-place to a DNA sequence by changing any Us to Ts
    """
    self.sequenceData = self.sequenceData.replace("U","T")

  def split(self, point = None):
    """
      Split this sequence into two halves and return them. The original sequence
      remains unmodified.

      :param point: defines the split point, if None then the centre is used
      :return: two Sequence objects -- one for each side
    """
    if point == None :
      point = len(self)/2

    r1 = FastqSequence(self.sequenceName + ".1",
                   self.sequenceData[:point])
    r2 = FastqSequence(self.sequenceName + ".2",
                   self.sequenceData[point:])

    return r1,r2

  def truncate(self, newLength):
    """
      Truncate this sequence in-place so it's only <newLength> nucleotides long.

      :param newLength: the length to truncate this sequence to.
    """
    return Sequence(self.sequenceName, self.sequenceData[:10])

  def clipThreePrime(self, seq, mm_score):
    """
      Clip a sequence from the 3' end of the sequence -- we assume the sequence
      to be clipped will always begin somewhere in this sequence, but may not be
      fully contained. If found, replaced with Ns.

      :param seq: sequence to be clipped
      :param mm_score: the number of matching bases needed to consider a hit,
                       mm_score = len(seq) would be 100% match
    """
    lim = mm_score - 1
    other_end = len(seq) - 1
    other_start = 0

    for i in range(len(self.sequenceData) - 1, lim - 1, -1) :
      self_end = i
      self_start = i - (len(seq) - 1)

      if self_start < 0 :
        self_start = 0
        other_start = other_end - self_end

      score = self.similarity(self_start, self_end, other_start, other_end, seq)
      if (score >= mm_score) :
        self.nsRight(len(seq) + (len(self) - i) - 1)
        break

  def clipAdaptor(self, adaptor):
    """
      Clip an adaptor sequence from this sequence. We assume it's in the 3' end.
      This is basically a convenience wrapper for clipThreePrime. It requires
      8 out of 10 of the first bases in the adaptor sequence to match for
      clipping to occur.

      :param adaptor: sequence to look for. We only use the first 10 bases; must
                      be a full Sequence object, not just a string.
    """
    missmatches = 2
    adaptor = adaptor.truncate(10)
    self.clipThreePrime(adaptor, len(adaptor) - missmatches)

  def containsAdaptor(self, adaptor):
    """
      Check whether this sequence contains adaptor contamination. If it exists,
      we assume it's in the 3' end. This function requires 8 out of 10 of the
      first bases in the adaptor sequence to match for an occurrence to be
      reported.

      :param adaptor: sequence to look for. We only use first 10 bases; must be
                      a full Sequence object, not just string.
      :return: True if there is an occurence of <adaptor>, False otherwise
    """
    origSeq = self.sequenceData
    self.clipAdaptor(adaptor)
    res = False
    if self.sequenceData != origSeq : res = True
    self.sequenceData = origSeq
    return res

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
    if len(mask) > len(self.sequenceData) : return False
    lim = len(mask)
    for i in range(0,lim):
      if mask[i] == "N" or mask[i] == "n" : continue
      if mask[i] != self.sequenceData[i] : return False
    return True


class FastaSequence(Sequence):
  """
    Data-structure for holding information about a Fasta-formatted sequence.
    This is basically just a wrapper around a regular Sequence object which
    provides the neccessary formatting for output.

    :param seqName:           the name of the sequence
    :param seqData:           the actual nucleotide sequence
    :param lineWidth:         width of sequence data lines for output. If None,
                              all sequence data will be output on a single line
    :param useMustableString:
  """

  def __init__(self, seqName, seqData = "", lineWidth = None,
               useMutableString = False):
    Sequence.__init__(self, seqName, seqData, useMutableString)
    self.lineWidth = lineWidth
  def __str__(self):
    """
      Get a string representation of this fasta sequence. If line width was
      specifified when the object was created, it is respected. Otherwise, all
      sequence data is printed to a single line.
    """
    if self.lineWidth == None :
      return ">" + self.sequenceName + "\n" + str(self.sequenceData)
    else :
      return self.formattedString()

  def formattedString(self):
    """
      Get a formatted version of the seq where no row of sequence data exceeds
      the line length for this object.

      :raise: SequenceError if this object has no line width specified.
    """
    if self.lineWidth == None :
      raise SequenceError("No line width for fasta formatted read specified")

    res = ">" + self.sequenceName + "\n"
    for i in range(0,len(self.sequenceData), self.lineWidth) :
      res += self.sequenceData[i:i+self.lineWidth]
      if i + self.lineWidth < len(self.sequenceData) : res += "\n"
    return res


class FastqSequence(Sequence):
  """
    Data structure fo holding informaton about a Fastq-formatted sequence.
    Fastq-formatted sequences differ from regular sequences by the inclusion of
    a quality score for each nucleotide in the sequence, encoded as a string.

    :param seqName:          A name describing the sequence. Can be any string.
    :param seqData:          The nucleotide sequence data. Can be DNA or RNA.
                             Note that there is no check to make sure the
                             sequence data is valid, that's the responsibility
                             of the caller.
    :param seqQual:          The quality string for this sequence -- must be
                             the same length as the nucleotide sequence.
    :param useMutableString: Store the sequence data as a mutable string, rather
                             than a regular python string. This should make
                             editing operations must faster, but it comes at the
                             expense of less flexibility (e.g. the object can
                             not be used as a hash key because it is mutable.)
    :raise SequenceError:    if the sequence data is not the same length as the
                             quality data.
  """

  def __init__(self, seqName, seqData=None,
               seqQual=None, useMutableString=False):
    """
      Constructor for FastqSequence class; see class level documentation for
      descriptions of parameters.
    """
    if len(seqData) != len(seqQual) :
      raise SequenceError("failed to create FastqSequence object -- length " +\
                          "of sequence data (" + str(len(seqData)) + ")" +\
                          "does not match length of quality string (" +\
                          str(len(seqQual)) + ")")

    Sequence.__init__(self, seqName, seqData, useMutableString)

    self.sequenceQual = seqQual

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
    return Sequence.__eq__(self, seq) and \
           self.sequenceQual == seq.sequenceQual

  def __ne__(self, read):
    """
      determine whether two fastqSequence objects are not equal. They are
      considered unequal if any of their sequence data (name, nuc. sequence)
      does not match, or if their quality data does not match.

      :param seq: the other sequence to compare against.
      :return: True if this sequence is not equal to seq, else False.
    """

    return Sequence.__ne__(self, read) or \
           self.sequenceQual != read.sequenceQual

  def truncate(self, size):
    """
      truncate this fastqSequence in-place so it is only <size> nucleotides long

      :param size: the number of nucleotides to truncate to.
    """
    self.trimRight(len(self) - size)

  def qualityToSolexa(self):
    """
      Convert the quality data for this fastqSequence in-place, from sanger to
      solexa format. Note that no checking is done to make sure the data was
      originally in sanger format; if it wasn't, the result will be junk.
    """
    newqual = ""
    for val in self.sequenceQual :
      newqual += chr(ord(val) + SANGER_SOLEXA_OFFSET)
    self.sequenceQual = newqual

  def trimRight(self, amount):
    """
      Trim this fastqSequence in-place by removing <amount> nucleotides from the
      3' end (right end).

      :param amount: the number of nucleotides to trim from the right-side of
                     this sequence.
    """
    self.sequenceData = self.sequenceData[:-amount]
    self.sequenceQual = self.sequenceQual[:-amount]

  def trimLeft(self, amount):
    """
      Trim this fastqSequence in-place by removing <amount> nucleotides from the
      5' end (left end).

      :param amount: the number of nucleotides to trim from the left-side of
                     this sequence.
    """
    self.sequenceData = self.sequenceData[amount:]
    self.sequenceQual = self.sequenceQual[amount:]

  def getRelativeQualityScore(self, i):
    val = self.sequenceQual[i]
    return (ord(val) - self.LOWSET_SCORE) / float (self.HIGHEST_SCORE -
                                                   self.LOWSET_SCORE)
  def reverseComplement(self):
    """
      Reverse complement this fastq sequence in-place.
    """
    Sequence.reverseComplement(self)
    self.sequenceQual = self.sequenceQual[::-1]

  def split(self, point = None):
    """
      Split this fastq sequence into two halves. The original sequence is left
      unaltered.

      :param point: the point (index) at which to split this sequence. If None
                    (the default), then we split in the middle.
      :return: two FastqSequence objects which correspond to the split of this
               sequence.
    """
    if point == None :
      point = len(self)/2

    r1 = FastqSequence(self.sequenceName + ".1",
              self.sequenceData[:point],
              self.sequenceQual[:point])
    r2 = FastqSequence(self.sequenceName + ".2",
              self.sequenceData[point:],
              self.sequenceQual[point:])
    return r1,r2

  def merge(self, other, forceMerge = False):
    """
      Merge two fastqSequences by concatenating their sequence data and their
      quality data (<self> first, then <other>); <self> and <other> must have
      the same sequence name. A new merged FastqSequence object is returned;
      <Self> and <other> are left unaltered.

      :param other: the other sequence to merge with self.
      :param forceMerge: force the merge to occur, even if sequences names don't
                         match. In this case, <self> takes precedence.
      :return: A new FastqSequence that represents the merging of <self> and
               <other>
      :raise: FastqSequenceError if the sequences names do not match, and the
              forceMerge parameter is not set.
    """
    if self.sequenceName != other.sequenceName and not forceMerge :
      raise FastqSequenceError("cannot merge " + self.sequenceName + " with " +\
                           other.sequenceName + " -- different sequence names")

    name = self.sequenceName
    seq = self.sequenceData + other.sequenceData
    qual = self.sequenceQual + other.sequenceQual

    return FastqSequence(name, seq, qual)

  def __str__(self):
    """
      Get a string representation of this FastQSequence.

      :return: String that represents this FastQSequence.
    """
    return "@" + self.sequenceName + "\n" + self.sequenceData +\
           "\n" + "+" + self.sequenceName + "\n" + self.sequenceQual


class SequenceUnitTests(unittest.TestCase):
  """
    Unit tests for sequence classes
  """

  def testClipadaptor(self):
    pass
    input =   Sequence("name", "ACTGCTAGCGATCGACT")
    adaptor = Sequence("adap",       "AGCGATAGACT")
    expect =  Sequence("name", "ACTGCTNNNNNNNNNNN")
    input.clipAdaptor(adaptor)
    got = input
    self.assertTrue(expect == got)

  def testNsLeft(self):
    input =   Sequence("name", "ACTGCTAGCGATCGACT")
    expect =  Sequence("name", "NNNNNTAGCGATCGACT")
    input.nsLeft(5)
    got = input
    self.assertTrue(expect == got)

  def testNsRight(self):
    input =   Sequence("name", "ACTGCTAGCGATCGACT")
    expect =  Sequence("name", "ACTGCTAGCGATNNNNN")
    input.nsRight(5)
    got = input
    self.assertTrue(expect == got)

  def testLengths(self):
    input =   Sequence("name", "ACTNCTANCGATNNACT")
    self.assertTrue(len(input) == 17)
    self.assertTrue(input.effectiveLength() == 13)

  def testMaskRegion(self):
    class TestRegion:
      def __init__(self, s, e):
        self.start = s
        self.end = e

    input =   Sequence("name", "ACTNCTANCGATNNACT")
    expect =  Sequence("name", "ANNNNTANCGATNNACT")
    out = input.copy()
    out.maskRegion(TestRegion(1,4))
    self.assertTrue(expect == out)

  def testCopyConstructor(self):
    one = Sequence("name", "ACTNCTANCGATNNACT")
    two = one.copy()
    self.assertTrue(one == two)
    one.sequenceName = "old"
    self.assertTrue(one != two)

  def testReverseComplement(self):
    input = Sequence("name", "ACTGCTAGCATGCGNN")
    expect = Sequence("name", "NNCGCATGCTAGCAGT")
    input.reverseComplement()
    self.assertTrue(input == expect)

  def testFormattedString(self):
    """
      test that string formatting works correctly for fasta sequences
    """
    r = FastaSequence("name", "ATCGATCGATCGATCTCGA", lineWidth=5)
    expect = ">name\n" +\
             "ATCGA\n" +\
             "TCGAT\n" +\
             "CGATC\n" +\
             "TCGA"
    got = r.formattedString()
    self.assertTrue(got == expect)

    # make sure this also works with a mutable underlying sequence
    r = FastaSequence("name", "ATCGATCGATCGATCTCGA", lineWidth=5,
                      useMutableString = True)
    expect = ">name\n" +\
             "ATCGA\n" +\
             "TCGAT\n" +\
             "CGATC\n" +\
             "TCGA"
    got = r.formattedString()
    self.assertTrue(got == expect)

  def testeq(self):
    """
      test the equality operator for fastQ sequences.
    """
    r1 = FastqSequence("s1","ACTGCT","BBBBBB")
    r2 = FastqSequence("s1","ACTGCT","BBBBBB")
    r3 = FastqSequence("s1","ACTGCT","BBBfBB")
    r4 = FastqSequence("s1","ACCGCT","BBBBBB")
    r5 = FastqSequence("s2","TCTGCT","fBBBBB")
    r6 = FastqSequence("s3","CCCCCC","fBBBBB")
    r7 = FastqSequence("s1","CCCCCC","fBBfBB")
    r7 = FastqSequence("s6","CCCCCC","fBBfBB")

    self.assertTrue((r1 == r2) == True)  # same name, same seq, same qual
    self.assertTrue((r1 == r3) == False) # same name, same seq, diff qual
    self.assertTrue((r1 == r4) == False) # same name, diff seq, same qual
    self.assertTrue((r3 == r4) == False) # same name, diff seq, diff qual
    self.assertTrue((r1 == r5) == False) # diff name, diff seq, diff qual
    self.assertTrue((r5 == r6) == False) # diff name, diff seq, same qual
    self.assertTrue((r6 == r7) == False) # diff name, same seq, diff qual
    self.assertTrue((r3 == r4) == False) # diff name, same seq, same qual

    self.assertTrue((r1 != r2) == False)  # same name, same seq, same qual
    self.assertTrue((r1 != r3) == True) # same name, same seq, diff qual
    self.assertTrue((r1 != r4) == True) # same name, diff seq, same qual
    self.assertTrue((r3 != r4) == True) # same name, diff seq, diff qual
    self.assertTrue((r1 != r5) == True) # diff name, diff seq, diff qual
    self.assertTrue((r5 != r6) == True) # diff name, diff seq, same qual
    self.assertTrue((r6 != r7) == True) # diff name, same seq, diff qual
    self.assertTrue((r3 != r4) == True) # diff name, same seq, same qual



if __name__ == "__main__":
    unittest.main()
