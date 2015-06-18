"""
  Date of Creation: 11th Dec 2014

  Description:   ...

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
import StringIO
import os
import os.path

# import for unit testing
import mock
import unittest

# pyokit imports
from pyokit.datastruct.genomeAlignment import GenomeAlignmentBlock
from pyokit.datastruct.genomeAlignment import GenomeAlignment
from pyokit.datastruct.sequence import Sequence
from pyokit.datastruct.sequence import UnknownSequence
from pyokit.io import maf


def genome_alignment_block_hash(b):
  """Hash a genome alignment block to unique ID: its location in the genome."""
  return b.chrom + "\t" + str(b.start) + "\t" + str(b.end)


###############################################################################
#    CLASSES FOR MANAGING ACCESS TO MULTIPLE GENOME ALIGNMENT FILES ON DISK   #
###############################################################################

def build_genome_alignment_from_directory(d_name, reference_species):
  """
  build a genome aligment by loading all files in a directory.

  Not recursive (i.e. subdirectories are not parsed). Will attempt to load all
  files regardless of extension. Expects all files to be MAF format genome
  alignment files.
  """
  blocks = []
  for fn in os.listdir(d_name):
    print fn
    pth = os.path.join(d_name, fn)
    print pth
    if os.path.isfile(pth):
      for b in genome_alignment_iterator(pth, reference_species):
        blocks.append(b)
  return GenomeAlignment(blocks)


###############################################################################
#                                ITERATORS                                    #
###############################################################################

def genome_alignment_iterator(fn, reference_species):
  """
  build an iterator for an MAF file of genome alignment blocks.

  :return an iterator that yields GenomeAlignment objects
  """
  kw_args = {"reference_species": reference_species}
  for e in maf.maf_iterator(fn, yield_class=GenomeAlignmentBlock,
                            yield_kw_args=kw_args):
    yield e


def load_jit_genome_alignment(alignment_fn, index_fn):
  """
  load a genome alignment from an index.

  :return: a just-in-time genome alignment, where access to any of the blocks
           causes them to be loaded from the alignment file
  """
  pass
  # factory = IndexFile(None, maf_iterator, genome_alignment_block_hash)


###############################################################################
#                                UNIT TESTS                                   #
###############################################################################

class TestGenomeAlignment(unittest.TestCase):
  def setUp(self):
    b1_hg19_seq = "atctccaagagggcataaaacac-tgagtaaacagctcttttatatgtgtttcctgga"
    b1_panTro_s = "atctccaagagggcataaaacac-tgagtaaacagctctt--atatgtgtttcctgga"
    b1_panTro_q = "99999999999999999999999-9999999999999999--9999999999999999"
    b1_tarSyr_s = "atctccaagagggctgaaaatgc-caaatga-----------tcacacgtttcctgga"
    b1_tarSyr_q = "79295966999999999999998-9999799-----------9999999999765775"
    b1_tupBel_s = "ttcaggaagggggcccaaaacgcttgagtggtcagctctta-ttttgcgtttactgga"
    b1_tupBel_q = "79648579699867994997775679665662767577569-6998745597677632"
    b1 = "a score=28680.000000\n" +\
         "s hg19.chr22             1711 57 + 51304566 " + b1_hg19_seq + "\n" +\
         "s panTro2.chrUn          1110 59 + 58616431 " + b1_panTro_s + "\n" +\
         "q panTro2.chrUn                             " + b1_panTro_q + "\n" +\
         "i panTro2.chrUn          C 0 C 0                          " + "\n" +\
         "s tarSyr1.scaffold_5923  2859 50 -     8928 " + b1_tarSyr_s + "\n" +\
         "q tarSyr1.scaffold_5923                     " + b1_tarSyr_q + "\n" +\
         "i tarSyr1.scaffold_5923  N 0 C 0                          " + "\n" +\
         "s tupBel1.scaffold_803   33686 61 +   85889 " + b1_tupBel_s + "\n" +\
         "q tupBel1.scaffold_803                      " + b1_tupBel_q + "\n" +\
         "i tupBel1.scaffold_803   I 1 C 0                          " + "\n" +\
         "e mm4.chr6            53310102 58 + 151104725 I"
    self.b1_hg19 = Sequence("hg19.chr22", b1_hg19_seq, 1711, 1768,
                            "+", 51302798)
    self.b1_panTro = Sequence("panTro2.chrUn", b1_panTro_s, 1110, 1169, "+",
                              58616431 - 1169,
                              {maf.QUALITY_META_KEY:b1_panTro_q,
                               maf.LEFT_STATUS_KEY:"C",
                               maf.LEFT_COUNT_KEY:0,
                               maf.RIGHT_STATUS_KEY:"C",
                               maf.RIGHT_COUNT_KEY:0})
    self.b1_tarSyr = Sequence("tarSyr1.scaffold_5923", b1_tarSyr_s,
                              8928 - 2859 - 50, 8928 - 2859, "-",
                              2859, {maf.QUALITY_META_KEY:b1_tarSyr_q,
                                     maf.LEFT_STATUS_KEY:"N",
                                     maf.LEFT_COUNT_KEY:0,
                                     maf.RIGHT_STATUS_KEY:"C",
                                     maf.RIGHT_COUNT_KEY:0})
    self.b1_mm4 = UnknownSequence("mm4.chr6", 53310102, 53310102 + 58, "+",
                                  151104725 - (53310102 + 58),
                                  {maf.EMPTY_ALIGNMENT_STATUS_KEY:"I"})

    b2_hg19_seq = "ccttcttttaattaattttgttaagg----gatttcctctagggccactgcacgtca"
    b2_panTro_s = "ccttcttttaattaattttgttatgg----gatttcgtctagggtcactgcacatca"
    b2_panTro_q = "99999999999999999999999999----999999099999999999999999999"
    b2_tarSyr_s = "tcttcttttaattaattttattgagggattgattccttattgggccactacacatta"
    b2_tarSyr_q = "999999899978999999999999999977989997998678865952859999899"
    b2_tupBel_s = "cct--gtttaaattactgtattg-gg----gatttcctatagggccgcttctcgtcc"
    b2_tupBel_q = "666--958759455555746366-68----656846556554745443677468565"
    b2 = "a score=31725.000000\n" +\
         "s hg19.chr22             1772 53 + 51304566 " + b2_hg19_seq + "\n" +\
         "s panTro2.chrUn          1169 53 + 58616431 " + b2_panTro_s + "\n" +\
         "q panTro2.chrUn                             " + b2_panTro_q + "\n" +\
         "i panTro2.chrUn          C 0 C 0                          " + "\n" +\
         "s tarSyr1.scaffold_5923  2909 124 -    8928 " + b2_tarSyr_s + "\n" +\
         "q tarSyr1.scaffold_5923                     " + b2_tarSyr_q + "\n" +\
         "i tarSyr1.scaffold_5923  C 0 N 0                          " + "\n" +\
         "s tupBel1.scaffold_803   33747 113 +  85889 " + b2_tupBel_s + "\n" +\
         "q tupBel1.scaffold_803                      " + b2_tupBel_q + "\n" +\
         "i tupBel1.scaffold_803 C 0 N 0              "
    self.maf1 = b1 + "\n\n" + b2
    self.b2_hg19 = Sequence("hg19.chr22", b2_hg19_seq, 1772, 1825,
                            "+", 51302741)
    self.b2_panTro = Sequence("panTro2.chrUn", b2_panTro_s, 1169, 1169 + 53,
                              "+", 58616431 - (1169 + 53),
                              {maf.QUALITY_META_KEY:b2_panTro_q,
                               maf.LEFT_STATUS_KEY:"C",
                               maf.LEFT_COUNT_KEY:0,
                               maf.RIGHT_STATUS_KEY:"C",
                               maf.RIGHT_COUNT_KEY:0})
    self.b2_tarSyr = Sequence("tarSyr1.scaffold_5923", b2_tarSyr_s,
                              8928 - 2909 - 124, 8928 - 2909, "-",
                              2909, {maf.QUALITY_META_KEY:b2_tarSyr_q,
                                     maf.LEFT_STATUS_KEY:"C",
                                     maf.LEFT_COUNT_KEY:0,
                                     maf.RIGHT_STATUS_KEY:"N",
                                     maf.RIGHT_COUNT_KEY:0})
    self.b1 = b1
    self.b2 = b2

  def test_genome_alignment_iterator(self):
    g_iter = genome_alignment_iterator(StringIO.StringIO(self.maf1), "hg19")
    blocks = [x for x in g_iter]

    # first test that everything for regular maf stuff works properly...
    self.assertEqual(len(blocks), 2)
    self.assertEqual(blocks[0]["hg19.chr22"], self.b1_hg19)
    self.assertEqual(blocks[0]["panTro2.chrUn"], self.b1_panTro)
    self.assertEqual(blocks[0]["tarSyr1.scaffold_5923"], self.b1_tarSyr)
    self.assertEqual(blocks[0]["mm4.chr6"], self.b1_mm4)
    self.assertEqual(blocks[1]["hg19.chr22"], self.b2_hg19)
    self.assertEqual(blocks[1]["panTro2.chrUn"], self.b2_panTro)
    self.assertEqual(blocks[1]["tarSyr1.scaffold_5923"], self.b2_tarSyr)

    self.assertEqual(maf.alignment_to_maf(blocks[0]).split(), self.b1.split())
    self.assertEqual(maf.alignment_to_maf(blocks[1]).split(), self.b2.split())
    out = StringIO.StringIO()
    maf.write_maf(genome_alignment_iterator(StringIO.StringIO(self.maf1),
                                            "hg19"),
                  out)
    self.assertEqual(self.maf1.split(), out.getvalue().split())

    # now we check to see that the extra stuff granted by genome alignments
    # works...
    self.assertEqual([b.start for b in blocks], [1711, 1772])
    self.assertEqual([b.end for b in blocks], [1711 + 57, 1772 + 53])
    self.assertEqual([b.chrom for b in blocks], ["chr22", "chr22"])
    self.assertEqual([b.reference_species for b in blocks], ["hg19", "hg19"])
    self.assertEqual([b.reference_sequence_name for b in blocks],
                     ["hg19.chr22", "hg19.chr22"])

  @mock.patch('os.listdir')
  @mock.patch('os.path.isfile')
  @mock.patch('__builtin__.open')
  def test_build_genome_alignment_from_directory(self, mock_open, mock_isfile,
                                                 mock_listdir):
    mock_listdir.return_value = ["one.maf", "two.maf", "some_sub_dir"]
    def open_side_effect(*args, **kwargs):
      if args[0] == os.path.join("the_dir", "one.maf"):
        return StringIO.StringIO(self.b1)
      elif args[0] == os.path.join("the_dir", "two.maf"):
        return StringIO.StringIO(self.b2)
      raise IOError("No such file")
    def isfile_side_effect(*args, **kwargs):
      if (args[0] == os.path.join("the_dir", "one.maf") or
          args[0] == os.path.join("the_dir", "two.maf")):
        return True
      return False
    mock_open.side_effect = open_side_effect
    mock_isfile.side_effect = isfile_side_effect

    ga = build_genome_alignment_from_directory("the_dir", "hg19")
    self.assertEqual(ga.num_blocks, 2)
    self.assertEqual(ga.get_blocks("chr22", 1711, 1720)[0]["hg19.chr22"],
                                   self.b1_hg19)
    self.assertEqual(ga.get_blocks("chr22", 1770, 1780)[0]["hg19.chr22"],
                                   self.b2_hg19)


###############################################################################
#              MAIN ENTRY POINT WHEN RUN AS STAND-ALONE MODULE                #
###############################################################################

if __name__ == "__main__":
    unittest.main()
