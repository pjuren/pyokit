#!/usr/bin/python

"""
  Date of Creation: 18th Dec 2014
  Description:      Data structures and functions for manipulating GO data

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


###############################################################################
#                          GENE ONTOLOGY CLASS                                #
###############################################################################

class GeneOntologyTerm(object):
  """
    Represents a gene ontology term, which may or may not have been scored
    for significance
  """

  def __init__(self, name, catagory=None):
    self.name = name
    self.catagory = catagory
    self.pvalue = None
    self.corrected_pvalue = None
