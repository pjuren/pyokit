#!/usr/bin/python

""" 
  Date of Creation: 24th November 2010     
                      
  Description:        Classes and functions for manipulating refseq 
                      transcriptome references

  Copyright (C) 2010  
  University of Southern California,
  Philip J. Uren,
  Andrew D. Smith
  
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
  
  --------------------
  
  Known Bugs:    None
  
  Revision 
  History:       None 
  
  
  TODO:          None 
"""

import random
from mapping.bedIterators import BEDIterator
from mapping.bed import toGenomicCoordinates

KEY_SEP = "___"

class TranscriptError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)

class Transcript :
  """
    @summary: a transcript is a collection of exons with a single refseq
              name
  """
  
  def __init__(self, exons):
    self.exons = exons
    self.start = min([exon.start for exon in exons])
    nms = [exon.name.split("_exon_")[0] for exon in exons]
    if len(set(nms)) != 1 : raise TranscriptError("exons must all be " +\
                                                  "from the same transcript: " + str(set(nms)))
    self.name = list(set(nms))[0]
    self.length = sum([exon.end - exon.start for exon in self.exons])
    chrms = [exon.chrom for exon in exons]
    if len(set(chrms)) != 1 : raise TranscriptError("exons must all be " +\
                                                  "from the same chromosome")
    self.chrom = list(set(chrms))[0]
    
  def generateRead(self, readLength):
    """
      @summary: generate a read in BEDElement from this transcript in 
                genomic coordinates
      @return: list of reads -- will be only 1 if the region covered 
               is only a single exon, otherwise it'll be a region for
               each exon covered
    """ 
    start = int(random.random() * (self.length - readLength))
    end = start + readLength
    return toGenomicCoordinates(start, end, self.exons)
  
  def transcriptID(self):
    return self.chrom + KEY_SEP + self.name

  def __len__(self):
    return self.length

def readReferenceAsTranscripts(reffn, verbose = False):
  transcripts = []
  ts = readTranscriptomeByGene(reffn, verbose = verbose)
  for k in ts :
    transcripts.append(Transcript(ts[k]))
  return transcripts
  

def makeReferenceKey(element):
  return element.chrom + KEY_SEP + element.name.split("_exon_")[0]

def getUniqueTranscriptIDs(elements):
  """
    @summary: given a list of exons taken from a transcriptome reference,
              get the list of chrom and refseq keys that uniquely ID
              each element
  """
  unique = set([element.chrom + KEY_SEP + element.name.split("_exon_")[0] 
                for element in elements])
  return list(unique)

def readTranscriptomeByGene(reffn, verbose):
  """
    @summary: read a transcriptome reference and group elements by gene
    @return: returns a dictionary where each item is a list of exons
             for a given gene. Entries are indexed by keys which are made
             by concatenating the chromosome and refseq separated by 
             KEY_SEP.
  """
  t = readTranscriptomeFlat(reffn, verbose = verbose)
  return splitReferenceByGene(t)
  
def splitReferenceByGene(elements, verbose = False):
  """
    @summary: take a list of elements from a transcritome reference 
              and split into lists for each refseq name, collected
              into a dictionary. Entries are indexed by keys which are made
              by concatenating the chromosome and refseq separated by 
              KEY_SEP.
  """
  transcripts = {}
  for element in elements :
    refseq = element.name.split("_exon_")[0]
    chrom = element.chrom
    key = chrom + KEY_SEP + refseq
    
    if not key in transcripts : transcripts[key] = []
    transcripts[key].append(element)
  return transcripts

def readTranscriptomeFlat(reffn, verbose = False):
  """
    @summary: read a transcriptome reference and return all elements
              as a flat list
    @return: list of all elements
  """
  return [element for element in BEDIterator(reffn, verbose = verbose)]
  
  