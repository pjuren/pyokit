#!/usr/bin/python

"""
Date of Creation: 23rd Feb 2014.

Description:      Register modules to run unit tests on.

Copyright (C) 2014
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
have_functioning_rpy2 = True
try:
  from rpy2.robjects import r
except ImportError, e:
  have_functioning_rpy2 = False

###############################################################################
#                           STATISTICS LIBRARIES                              #
###############################################################################
if have_functioning_rpy2:
  from pyokit.statistics.multipleHypothesisTesting import TestMHT
  from pyokit.statistics.fisher import FisherTests
from pyokit.statistics.online import TestOnlineStats
from pyokit.statistics.probability import ProbabilityTests
from pyokit.statistics.beta import BetaDistTests

###############################################################################
#                              I/O LIBRARIES                                  #
###############################################################################
from pyokit.io.bedIterators import BEDIteratorUnitTests
from pyokit.io.fastqIterators import FastQUintTests
from pyokit.io.fastaIterators import TestFastaIterators
from pyokit.io.repeatmaskerAlignments import TestAlignmentIterators
from pyokit.io.indexedFile import TestIndexedFile
from pyokit.io.repeatmaskerAnnotations import TestRepMaskerIterators
from pyokit.io.maf import TestMAF
from pyokit.io.genomeAlignment import TestGenomeAlignment

from pyokit.io.wigDir import *
from pyokit.io.wigFile import *
from pyokit.io.indexedWig import *
from pyokit.io.wigIterators import WigIteratorUnitTests

###############################################################################
#                              DATA STRUCTURES                                #
###############################################################################
from pyokit.datastruct.multipleAlignment import TestAlignments
from pyokit.datastruct.retrotransposon import TestRetrotransposon
from pyokit.datastruct.read import NGSReadUnitTests
from pyokit.datastruct.genomeAlignment import TestGenomeAlignmentDS
from pyokit.datastruct.intervalTree import *
from pyokit.datastruct.genomicInterval import TestGenomicInterval
from pyokit.datastruct.sequence import *

###############################################################################
#                             TESTING AND UTILS                               #
###############################################################################
from pyokit.testing.dummyfiles import *
from pyokit.util.progressIndicator import TestProgressIndicator
from pyokit.util.meta import TestMeta

###############################################################################
#                                  SCRIPTS                                    #
###############################################################################
if have_functioning_rpy2:
  from pyokit.scripts.fdr import TestFDR
from pyokit.scripts.conservationProfile import TestConservationProfileIndvFiles
from pyokit.scripts.conservationProfile import TestConservationProfileDirectory
from pyokit.scripts.index import TestIndex
from pyokit.scripts.join import TestJoin
from pyokit.scripts.fastqNucDist import TestNucDist
from pyokit.scripts.remDupsBED import TestRemDupsBed
from pyokit.scripts.convertJunctionReads import ConvertJunctionsUnitTests
from pyokit.scripts.regionCollapse import TestCollapseRegions
from pyokit.scripts.overlapProfile import TestOverlapProfile
from pyokit.scripts.pyokitMain import PyokitSmokeTests


if __name__ == "__main__":
  head = " ------------------------  DATA STRUCTURES  --------------------- \n"
  sys.stderr.write(head)
  sys.stderr.write("registered tests in \n")
  sys.stderr.write("  " + str(NGSReadUnitTests) + "\n")
  sys.stderr.write("  " + str(TestAlignments) + "\n")
  sys.stderr.write("  " + str(TestRetrotransposon) + "\n")
  sys.stderr.write("  " + str(TestGenomicInterval) + "\n")
  sys.stderr.write("  " + str(TestGenomeAlignmentDS) + "\n\n")

  head = " ------------------------        IO         --------------------- \n"
  sys.stderr.write(head)
  sys.stderr.write("registered tests in \n")
  sys.stderr.write("  " + str(TestAlignmentIterators) + "\n")
  sys.stderr.write("  " + str(TestIndexedFile) + "\n")
  sys.stderr.write("  " + str(BEDIteratorUnitTests) + "\n")
  sys.stderr.write("  " + str(FastQUintTests) + "\n")
  sys.stderr.write("  " + str(TestFastaIterators) + "\n")
  sys.stderr.write("  " + str(WigIteratorUnitTests) + "\n")
  sys.stderr.write("  " + str(TestRepMaskerIterators) + "\n")
  sys.stderr.write("  " + str(TestMAF) + "\n")
  sys.stderr.write("  " + str(TestGenomeAlignment) + "\n\n")

  head = " ------------------------     STATISTICS    --------------------- \n"
  sys.stderr.write(head)
  sys.stderr.write("registered tests in \n")
  sys.stderr.write("  " + str(ProbabilityTests) + "\n")
  sys.stderr.write("  " + str(TestOnlineStats) + "\n")
  sys.stderr.write("  " + str(BetaDistTests) + "\n")
  if have_functioning_rpy2:
    sys.stderr.write("  " + str(TestMHT) + "\n")
    sys.stderr.write("  " + str(FisherTests) + "\n")
  sys.stderr.write("\n")

  head = " ------------------------  TESTING AND UTIL  -------------------- \n"
  sys.stderr.write(head)
  sys.stderr.write("registered tests in \n")
  sys.stderr.write("  " + str(TestProgressIndicator) + "\n")
  sys.stderr.write("  " + str(TestMeta) + "\n\n")

  head = " -----------------------------  SCRIPTS  ------------------------ \n"
  sys.stderr.write(head)
  sys.stderr.write("registered tests in \n")
  if have_functioning_rpy2:
    sys.stderr.write("  " + str(TestFDR) + "\n")
  sys.stderr.write("  " + str(TestIndex) + "\n")
  sys.stderr.write("  " + str(TestJoin) + "\n")
  sys.stderr.write("  " + str(TestConservationProfileDirectory) + "\n")
  sys.stderr.write("  " + str(TestCollapseRegions) + "\n")
  sys.stderr.write("  " + str(TestNucDist) + "\n")
  sys.stderr.write("  " + str(TestRemDupsBed) + "\n")
  sys.stderr.write("  " + str(ConvertJunctionsUnitTests) + "\n")
  sys.stderr.write("  " + str(TestConservationProfileIndvFiles) + "\n")
  sys.stderr.write("  " + str(TestOverlapProfile) + "\n")
  sys.stderr.write("  " + str(PyokitSmokeTests) + "\n\n")

  sys.stderr.write("\n\n RUNNING TESTS \n\n")
  unittest.main()
