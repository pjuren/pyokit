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
import os
import sys

# pyokit imports
from pyokit.datastruct import retrotransposon
from pyokit.datastruct.multipleAlignment import JustInTimePairwiseAlignment
from pyokit.datastruct import multipleAlignment
from pyokit.datastruct.genomicInterval import GenomicInterval
from pyokit.io.alignmentIterators import repeat_masker_alignment_iterator
from pyokit.io.indexedFile import IndexedFile
from pyokit.io.indexedFile import IndexError
from pyokit.util.progressIndicator import ProgressIndicator


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

  # try to get an idea of how much data we have...
  if verbose :
    try :
      total = os.path.getsize(strm.name)
      pind = ProgressIndicator(totalToDo=total, messagePrefix="completed",
                               messageSuffix="of processing " + strm.name)
    except AttributeError as e:
      sys.stderr.write(str(e))
      sys.stderr.write("completed [unknown] of processing index")
      verbose = False

  if header:
    # chomp first 2 lines
    next(strm)
    next(strm)

  for line in strm:
    if verbose :
      pind.done = strm.tell()
      pind.showProgress()

    line = line.strip()
    if line == "":
      continue
    rto = retrotransposon.from_repeat_masker_string(line)
    if alignment_index != None:
      rto.pairwise_alignment =\
          JustInTimePairwiseAlignment(alignment_index, rto.uniq_id)
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
    an1 = "283 26.37 4.21 0.00 chr1 15 67 (266) C B A (119) 141 85 5"
    an2 = "318 22.50 3.61 0.00 chr1 15266 15323 (249235276) C " +\
          "MIR3 SINE/MIR (65) 143 84 10"
    an3 = "18 23.18 0.00 1.96 chr1 15798 15830 (249234772) " +\
          "(TGCTCC)n Simple_repeat 1 32 (0) 15"
    an4 = "487 20.75 0.93 0.93 chr1 158389 158409 (249092126) C " +\
          "Charlie29b DNA/hAT-Charlie (532) 662 641 231"
    an5 = "318  23.0  3.7  0.0  chr1        15265   15355 (249235266) C  " +\
          "MIR3           SINE/MIR             (119)  143     49      300"
    an6 = "18  23.2  0.0  2.0  chr1        15798   15849 (249234772) +  " +\
          "(TGCTCC)n      Simple_repeat            1   51    (0)      319"
    self.ann = (head_l1 + "\n" + head_l2 + "\n\n" + an1 + "\n" + an2 + "\n" +
                an3 + "\n" + an4 + "\n" + an5 + "\n" + an6)
    self.indv_an = [an1, an2, an3, an4, an5, an6]

    # set up a repeat-masker alignment file...
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
               "Gap_init rate = 0.03 (2 / 79), avg. gap size = 1.50 (3 / 2)"

    alig_3_header = "18 23.18 0.00 1.96 chr1 15798 15830 (249234772) " +\
                    "(TGCTCC)n#Simple_repeat 1 32 (0) m_b1s252i0 15"
    alig_3 = "  chr1          15798 GCTGCTTCTCCAGCTTTCGCTCCTTCATGCT 15828\n" +\
             "                         v  v    v   iii      v - v        \n" +\
             "  (TGCTCC)n#Sim     1 GCTCCTGCTCCTGCTCCTGCTCCTGC-TCCT 30   \n" +\
             "                                                           \n" +\
             "  chr1          15829 GC 15830                             \n" +\
             "                                                           \n" +\
             "  (TGCTCC)n#Sim    31 GC 32                                  "
    alig_3_m = "Matrix = Unknown                                   \n" +\
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

  def test_basic_iterator(self):
    elems = [x for x in repeat_masker_iterator(StringIO.StringIO(self.ann))]
    self.assertEqual(len(elems), 6)
    for i in range(0, len(elems)):
      an = retrotransposon.from_repeat_masker_string(self.indv_an[i])
      self.assertEqual(elems[i], an)

    # alignments are not avaialble, so liftover should work only on coords;
    # just check one to make sure its working
    # elem[0]: 15, 67 -> 85, 141 - (53 to 57; gap_length = 14)
    self.assertEqual(elems[0].liftover(GenomicInterval("chr1", 10, 100)),
                     [GenomicInterval("A#B", 128, 142, strand='+'),
                      GenomicInterval("A#B", 113, 127, strand='+'),
                      GenomicInterval("A#B", 98, 112, strand='+'),
                      GenomicInterval("A#B", 86, 97, strand='+')])

  def test_iterator_with_alignment_index(self):
    def extract_UID(rm_alignment):
      return rm_alignment.meta[multipleAlignment.RM_ID_KEY]

    s_io = StringIO.StringIO(self.rm_rc_1_input)
    index = IndexedFile(s_io, repeat_masker_alignment_iterator, extract_UID)

    elems = [x for x in repeat_masker_iterator(StringIO.StringIO(self.ann),
                                               alignment_index=index)]
    self.assertEqual(len(elems), 6)
    for i in range(0, len(elems)):
      an = retrotransposon.from_repeat_masker_string(self.indv_an[i])
      self.assertEqual(elems[i], an)

    # alignments were provided, liftover should be using them; test one
    # to make sure they were matched up properly
    r = elems[0].liftover(GenomicInterval("chr1", 10, 100))
    self.assertEqual(r, [(132, 142), (120, 131), (88, 118), (85, 87)])
    # also test one of the ones that had no alignment; here we expect failure
    self.assertRaises(IndexError, elems[4].liftover,
                      GenomicInterval("chr1", 15200, 15400))



###############################################################################
#                 ENTRY POINT IF RUN AS STAND-ALONE MODULE                    #
###############################################################################

if __name__ == '__main__':
    unittest.main()
