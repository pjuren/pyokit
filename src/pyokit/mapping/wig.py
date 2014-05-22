#!/usr/bin/python

""" 
  Date of Creation: 11th November 2010     
                       
  Description:   Classes and functions for manipulating wiggle format
                 data

  Copyright (C) 2010  
  University of Southern California,
  Philip J. Uren,
  Andrew D. Smith
  
  Authors: Philip J. Uren, Andrew D. Smith
  
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
  History:       18th Oct 2011 -- Philip J. Uren
                   * added < operator
  
  TODO:          * Class and method headers are missing
                 * Allow choice of ordering by end 
"""

import sys, unittest
import weakref

class WigError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)

class WigElement(object):
  # define the WigElements slots -- saves us lots of memory
  # waste with the dynamic dictionary for names!
  __slots__ = ['chrom', 'start', 'end', 'score']
  
  def __init__(self, chrom, start, end, score):
    """
      @summary: constructor for WigElement objects
    """ 
    self.chrom = chrom
    self.start = start
    self.end = end
    self.score = score
  
  def __str__(self):
    """
      @summary: get string representation of WigElement object
    """
    if float(int(self.score)) == self.score : sc = str(int(self.score))
    else : sc = str(self.score) 
    return self.chrom + "\t" + str(self.start) + "\t" + str(self.end) +\
           "\t" + sc
           
  def before(self, w2) :
    """
      @summary: return true if self is before w2
      @deprecated: use < operator 
    """
    if w2 == None : return False
    if self.chrom < w2.chrom or (self.chrom == w2.chrom and \
                                 self.start < w2.start) : return True
    return False
  
  def __lt__(self, rhs) :
    """
      @summary: less than operator for WigElements
    """
    return self.before(rhs)

def wigElementFromString(s, scoreType=int):
  parts = s.split("\t")
  if (len(parts) < 4) :
    raise WigError("failed to parse " + s + " as wig element, too few fields")
  return WigElement(parts[0].strip(), int(parts[1]), int(parts[2]), 
                    scoreType(parts[3]))


if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])