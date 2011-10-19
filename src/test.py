#!/usr/bin/python
import unittest

########### MAPPING LIBRARIES ############
from mapping.bed import *          
from mapping.bedIterators import * 

from mapping.wig import *
from mapping.wigDir import *
from mapping.wigFile import *
from mapping.indexedWig import *
from mapping.wigIterators import *

######### SEQUENCING LIBRARIES ###########
from sequencing.fasta import *
from sequencing.fastaread import *
from sequencing.fastq import *
from sequencing.fastqread import *
from sequencing.fastread import *
from sequencing.maf import *

############ DATA STRUCTURES #############
from datastruct.intervalTree import *

########### TESTING AND UTILS ############
from testing.dummyfiles import *
from util.progressIndicator import *


if __name__ == "__main__":
  unittest.main()

