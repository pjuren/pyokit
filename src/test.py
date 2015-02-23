#!/usr/bin/python

"""
  Date of Creation: 23rd Feb 2014
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

###############################################################################
#                           STATISTICS LIBRARIES                              #
###############################################################################
from pyokit.statistics.probability import ProbabilityTests

###############################################################################
#                              I/O LIBRARIES                                  #
###############################################################################
from pyokit.io.bedIterators import BEDIteratorUnitTests
from pyokit.io.fastqIterators import FastQUintTests
from pyokit.io.fastaIterators import TestFastaIterators
from pyokit.io.alignmentIterators import TestAlignmentIterators
from pyokit.io.indexedFile import TestIndexedFile

from pyokit.io.wigDir import *
from pyokit.io.wigFile import *
from pyokit.io.indexedWig import *
from pyokit.io.wigIterators import WigIteratorUnitTests

###############################################################################
#                              DATA STRUCTURES                                #
###############################################################################
from pyokit.datastruct.multipleAlignment import TestAlignments
from pyokit.datastruct.retrotransposon import TestRetrotransposon
from pyokit.datastruct.intervalTree import *
from pyokit.datastruct.genomicInterval import *
from pyokit.datastruct.sequence import *
from pyokit.datastruct.maf import *

###############################################################################
#                             TESTING AND UTILS                               #
###############################################################################
from pyokit.testing.dummyfiles import *
from pyokit.util.progressIndicator import *


if __name__ == "__main__":
  sys.stderr.write("registered tests in " + str(TestAlignmentIterators) + "\n")
  sys.stderr.write("registered tests in " + str(TestIndexedFile) + "\n")
  sys.stderr.write("registered tests in " + str(BEDIteratorUnitTests) + "\n")
  sys.stderr.write("registered tests in " + str(FastQUintTests) + "\n")
  sys.stderr.write("registered tests in " + str(TestFastaIterators) + "\n")
  sys.stderr.write("registered tests in " + str(WigIteratorUnitTests) + "\n")
  sys.stderr.write("registered tests in " + str(TestRetrotransposon) + "\n")
  sys.stderr.write("registered tests in " + str(TestAlignments) + "\n")
  sys.stderr.write("registered tests in " + str(ProbabilityTests) + "\n")
  unittest.main()
