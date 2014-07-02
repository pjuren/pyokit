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
  DNA_COMPLEMENTS = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N",
                     "a":"t", "t":"a", "c":"g", "g":"c", "n":"n"}
  RNA_COMPLEMENTS = {"A":"U", "U":"A", "C":"G", "G":"C", "N":"N",
                     "a":"u", "u":"a", "c":"g", "g":"c", "n":"n"}

  def __init__(self, seqName, seqData, useMutableString = False):
    """
      @summary: Constructor
      @param seqName:
      @param seqData:
      @param useMutableString:
    """
    self.sequenceName = seqName
    if useMutableString :
      self.sequenceData = MutableString(seqData)
    else :
      self.sequenceData = seqData
    self.mutableString = useMutableString

  def copy(self):
    """
      @summary: Copy constructor
    """
    return Sequence(self.sequenceName, self.sequenceData, self.mutableString)

  def percentNuc(self, nuc):
    """
      @summary: return the percentage of the sequence which is equal
              to the passed nucleotide
      @param nuc: count number of times this nuc appears
      @return: percentage (float)
    """
    count = reduce(lambda x,y: x+1 if y==nuc else x, self.sequenceData, 0)
    return count / float(len(self.sequenceData))

  def similarity(self, self_start, self_end, other_start, other_end, other):
    """
      DESCPT: Count of the number of places where this[start,end] is equal to
              other[o_start, o_end]
    """
    assert(self_end - self_start == other_end - other_start)
    count = 0
    for i in range(0, self_end - self_start + 1) :
      if self.sequenceData[self_start + i] == other.sequenceData[other_start + i] :
        count += 1
    return count

  def reverseComplement(self):
    """
      @summary: reverse complement this sequence
    """
    isRNA = self.isRNA()
    tmp = ""
    for n in self.sequenceData :
      if isRNA : tmp += Sequence.RNA_COMPLEMENTS[n]
      else : tmp += Sequence.DNA_COMPLEMENTS[n]
    self.sequenceData = tmp[::-1]

  def __len__(self):
    """
      DESCPT: length of the read, defined as the length of its sequence data
    """
    return len(self.sequenceData)

  def effectiveLength(self):
    """
      DESCPT: disregarding N's, how long is this sequence?
    """
    return len([nuc for nuc in self.sequenceData
                    if nuc != "N" and nuc != "n"])

  def __eq__(self, read):
    """
      DESCRP: return true if this read is equal to passed parameter
    """
    if read == None : return False
    return  self.sequenceData == read.sequenceData and\
            self.sequenceName == read.sequenceName

  def __ne__(self, read):
    """
      DECPRP: return true if this read is not equal to the passed parameter
    """
    if read == None : return True
    return  self.sequenceData != read.sequenceData or\
            self.sequenceName != read.sequenceName

  def nsLeft(self, amount):
    """
      DESCPT: replace leftmost <amount> bases by Ns
    """
    self.sequenceData = (amount * "N") + self.sequenceData[amount:]

  def nsRight(self, amount):
    """
      DESCPT: replace rightmost <amount> bases by Ns
    """
    self.sequenceData = self.sequenceData[:-amount] + (amount * "N")

  def maskRegion(self, region):
    """
      DESCPT: Replace nucleotides in this read in the regions
              given by Ns
      PARAMS: region -- any object with .start and .end attributes
                        co-ords are zero based and inclusive of both
                        end points
      RAISES: SequenceError -- if region specifies nucleotides not present in
                               this read
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
      DESCPT: Mask the given regions in this read
      PARAMS: verbose -- print status messages to stderr if True
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
      @summary: return True if this sequence contains only DNA nucleotides
      @return: True if contains only DNA nucleotides, False otherwise
    """
    for nuc in self.sequenceData :
      if not nuc in "ACGTacgtn" : return False
    return True

  def isRNA(self):
    """
      @summary: return True if this sequence contains only RNA nucleotides
      @return: True if contains only RNA nucleotides, False otherwise
    """
    for nuc in self.sequenceData :
      if not nuc in "ACGUacgun" : return False
    return True

  def toRNA(self):
    """
      @summary: convert to RNA sequence by changing any Ts to Us
    """
    self.sequenceData = self.sequenceData.replace("T","U")

  def toDNA(self):
    """
      @summary: convert to DNA sequence by changing any Us to Ts
    """
    self.sequenceData = self.sequenceData.replace("U","T")

  def split(self, point = None):
    """
      DESCPT: Split this read into two halves
      PARAMS: point -- defines the split point, if None then the centre is used
      RETURN: two Sequence objects -- one for each side
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
      DESCRP: truncate the read so it is only <newLength> nucleotides long
    """
    return Sequence(self.sequenceName, self.sequenceData[:10])

  def clipThreePrime(self, seq, mm_score):
    """
      DESCPT: Clip a sequence from the 3' end of the read -- we assume sequence to
              be clipped will always begin somewhere in this sequence, but may
              not be fully contained. If found, replaced with Ns
      PARAMS: seq -- sequence to be clipped
              mm_score -- the number of matching bases needed to consider a hit,
                          mm_score = len(adaptor) would be 100% match
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
      DESCPT: Clip adaptor sequence from this read. We assume it's in the
              3' end
      PARAMS: adaptor -- sequence to look for. We only use first 10 bases
                         must be a full Sequence object, not just string
    """
    missmatches = 2
    adaptor = adaptor.truncate(10)
    self.clipThreePrime(adaptor, len(adaptor) - missmatches)

  def containsAdaptor(self, adaptor):
    """
      DESCPT: Does this sequence contain adaptor contamination?
              We assume adaptor is in 3' end
      PARAMS: adaptor -- sequence to look for. must be a full Sequence
                         object, not just string
      RETURN: bool -- true if there is an occurence of <adaptor>,
              false otherwise
    """
    origSeq = self.sequenceData
    self.clipAdaptor(adaptor)
    res = False
    if self.sequenceData != origSeq : res = True
    self.sequenceData = origSeq
    return res

  def isPolyA(self):
    """
      DESCRP: Is this sequence a polyA? based on having > 90% A in read
    """
    return self.percentNuc("A") >= 0.9

  def isPolyT(self):
    """
      DESCRP: Is this sequence a polyT? based on having > 90% A in read
    """
    return self.percentNuc("T") >= 0.9

  def isLowQuality(self):
    """
      DESCRP: Is this sequence low quality? based on having > 10% N in read
    """
    return self.percentNuc("N") >= 0.1

  def maskMatch(self, mask):
    """
      DESCPT: Determine whether this sequence matches the given mask.
              Ns in the mask are considered to match anything in the
              sequence -- all other chars must match exactly
      PARAMS: mask -- string to match against
      RETURN: True if the mask matches at all places, otherwise false
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
  def __init__(self, seqName, seqData=None,
               seqQual=None, useMutableString=False):
    Sequence.__init__(self, seqName, seqData, useMutableString)
    self.sequenceQual = seqQual

    # for quality scores
    self.LOWSET_SCORE = 64
    self.HIGHEST_SCORE = 104

  def __eq__(self, read):
    return Sequence.__eq__(self, read) and \
           self.sequenceQual == read.sequenceQual

  def __ne__(self, read):
    return Sequence.__ne__(self, read) or \
           self.sequenceQual != read.sequenceQual

  def truncate(self, size):
    self.trimRight(len(self) - size)

  def qualityToSolexa(self):
    """
      @summary: convert quality data from sanger to solexa format. Note that
                no checking is done to make sure the data was originally in
                sanger format; if it wasn't, the result will be junk
    """
    newqual = ""
    for val in self.sequenceQual :
      newqual += chr(ord(val) + SANGER_SOLEXA_OFFSET)
    self.sequenceQual = newqual

  def trimRight(self, amount):
    """
      @summary: trim the read by removing <amount> nucleotides
                from the 3' end (right end)
    """
    self.sequenceData = self.sequenceData[:-amount]
    self.sequenceQual = self.sequenceQual[:-amount]

  def trimLeft(self, amount):
    self.sequenceData = self.sequenceData[amount:]
    self.sequenceQual = self.sequenceQual[amount:]

  def getRelativeQualityScore(self, i):
    val = self.sequenceQual[i]
    return (ord(val) - self.LOWSET_SCORE) / float (self.HIGHEST_SCORE -
                                                   self.LOWSET_SCORE)
  def reverseComplement(self):
    """
      Reverse complement this fastq sequence.
    """
    Sequence.reverseComplement(self)
    self.sequenceQual = self.sequenceQual[::-1]

  def split(self, point = None):
    """ returns two Read objects which correspond to the split of this read """
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
    if self.sequenceName != other.sequenceName and not forceMerge :
      raise FastqSequenceError("cannot merge " + self.sequenceName + " with " +\
                           other.sequenceName + " -- different sequence names")

    name = self.sequenceName
    seq = self.sequenceData + other.sequenceData
    qual = self.sequenceQual + other.sequenceQual

    return FastqSequence(name, seq, qual)

  def __str__(self):
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
    r = FastaSequence("name", "ATCGATCGATCGATCTCGA")
    expect = ">name\n" +\
             "ATCGA\n" +\
             "TCGAT\n" +\
             "CGATC\n" +\
             "TCGA"
    got = r.formattedString(width=5)
    self.assertTrue(got == expect)

    # make sure this also works with a mutable underlying sequence
    r = FastaSequence("name", "ATCGATCGATCGATCTCGA", useMutableString = True)
    expect = ">name\n" +\
             "ATCGA\n" +\
             "TCGAT\n" +\
             "CGATC\n" +\
             "TCGA"
    got = r.formattedString(width=5)
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
