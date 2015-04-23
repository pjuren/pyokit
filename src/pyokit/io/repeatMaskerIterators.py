#!/usr/bin/python

"""
  Date of Creation: 11th Dec 2014
  Description:      Iterators for processing RepeatMasker files

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

# Standard imports
import unittest
import StringIO

# pyokit imports
from pyokit.datastruct import retrotransposon
from pyokit.datastruct.multipleAlignment import JustInTimePairwiseAlignment


def repeat_masker_iterator(fh, alignment_index=None,
                           header=True, verbose=False):
  """
  Iterator for repeatmasker coordinate annotation files. These files describe
  the location of repeat occurrences. There is (optionally) a two-line header
  with the names of the fields (ignored by the iterator, if present). Each line
  is a record of an occurrence. The description of fields for each line is
  given in from_repeat_masker_string.

  :param fh:              stream-like object, or string filename, to load the
                          annotations from
  :param alignment_index: an IndexedFile for full alignments; keys should be
                          repeat-masker IDs
  :param header:          if True, expect and discard the two-line header;
                          otherwise we will expect there is no header
  :param verbose:         if True, output additional status messages about
                          progress to stderr.
  """

  strm = fh
  if type(fh).__name__ == "str":
    strm = open(fh)

  if header:
    # chomp first 2 lines
    next(strm)
    next(strm)

  for line in strm:
    line = line.strip()
    if line == "":
      continue
    rto = retrotransposon.from_repeat_masker_string(line)
    if alignment_index != None:
      rto.pairwise_alignment =\
          JustInTimePairwiseAlignment(alignment_index, rto.meta[RM_ID_KEY])
    yield rto


###############################################################################
#                                 UNIT TESTS                                  #
###############################################################################

class TestRepMaskerIterators(unittest.TestCase):

  def setUp(self):
    """
      Set up some data to use in tests; create a repeatmasker annotation file
      and a corresponding full alignment file.
    """

    # annotations are space separated; this is just verbatim copied from the
    # first few lines in the hg19 repeatmasker annotation
    head_l1 = "SW  perc perc perc  query      position in query           " +\
              "matching       repeat              position in  repeat"
    head_l2 = "score  div. del. ins.  sequence    begin     end    (left) " +\
              "   repeat         class/family         begin  end (left)   ID"
    an1 = "463   1.3  0.6  1.7  chr1        10001   10468 (249240153) +  " +\
          "(TAACCC)n      Simple_repeat            1  463    (0)      1"
    an2 = "3612  11.4 21.5  1.3  chr1        10469   11447 (249239174) C  " +\
          "TAR1           Satellite/telo       (399) 1712    483      2"
    an3 = "484  25.1 13.2  0.0  chr1        11505   11675 (249238946) C  " +\
          "L1MC5a         LINE/L1             (2382) 5648   5452      3"
    an4 = "239  29.4  1.9  1.0  chr1        11678   11780 (249238841) C" +\
          "  MER5B          DNA/hAT-Charlie       (74)  104      1      4"
    an5 = "318  23.0  3.7  0.0  chr1        15265   15355 (249235266) C  " +\
          "MIR3           SINE/MIR             (119)  143     49      5"
    an6 = "18  23.2  0.0  2.0  chr1        15798   15849 (249234772) +  " +\
          "(TGCTCC)n      Simple_repeat            1   51    (0)      6"
    self.ann = (head_l1 + "\n" + head_l2 + "\n\n" + an1 + "\n" + an2 + "\n" +
                an3 + "\n" + an4 + "\n" + an5 + "\n" + an6)
    self.indv_an = [an1, an2, an3, an4, an5, an6]

  def test_basic_iterator(self):
    elems = [x for x in repeat_masker_iterator(StringIO.StringIO(self.ann))]
    self.assertEqual(len(elems), 6)
    for i in range(0, len(elems)):
      an = retrotransposon.from_repeat_masker_string(self.indv_an[i])
      self.assertEqual(elems[i], an)

  def test_iterator_with_index(self):
    pass

###############################################################################
#                 ENTRY POINT IF RUN AS STAND-ALONE MODULE                    #
###############################################################################

if __name__ == '__main__':
    unittest.main()
