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
    Represents a gene ontology term.
  """

  def __init__(self, name, identified=None, catagory=None):
    """
    :param name:        the term name, can contain spaces;
                        e.g.: intracellular transport
    :param identifier:  the term identifier; e.g: GO:0046907
    :param catagory:    the database or catagory the term belongs to;
                        e.g.: GOTERM_BP_FAT
    """
    self.name = name
    self.catagory = catagory


class GeneOntologyEnrichmentResult(GeneOntologyTerm):
  """
  Represents the result of a gene ontology enrichment calculation for a
  single GO term.
  """

  def __init__(self, name, pvalue, identifier=None, catagory=None):
    """

    """
    super(GeneOntologyEnrichmentResult, self).__init__(name,
                                                       identifier, catagory)
    self.pvalue = pvalue
