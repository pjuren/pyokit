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

############################## I/O LIBRARIES ###################################
from pyokit.io.bedIterators import *
from pyokit.io.fastaIterators import *
from pyokit.io.fastqIterators import *

from pyokit.io.wigDir import *
from pyokit.io.wigFile import *
from pyokit.io.indexedWig import *
from pyokit.io.wigIterators import *

######### SEQUENCING LIBRARIES ###########
from pyokit.sequencing.maf import *

############ DATA STRUCTURES #############
from pyokit.datastruct.intervalTree import *
from pyokit.datastruct.genomicInterval import *
from pyokit.datastruct.sequence import *

########### TESTING AND UTILS ############
from pyokit.testing.dummyfiles import *
from pyokit.util.progressIndicator import *


if __name__ == "__main__":
  unittest.main()
