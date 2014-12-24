#!/usr/bin/python

"""
  Date of Creation: 19th Dec 2014
  Description:      Functions for loading DAVID gene ontology results

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

# pyokit imports
from pyokit.datastruct.geneOntology import GeneOntologyEnrichmentResult

###############################################################################
#                           MODULE-LEVEL CONSTANTS                            #
###############################################################################

NUM_FIELDS_IN_DAVID_RECORD = 13
PVAL_FIELD_NUM = 11


###############################################################################
#                                  ITERATORS                                  #
###############################################################################

def david_results_iterator(fn, verbose=False):
  """
  Iterate over a DAVID result set and yeild GeneOntologyTerm objects
  representing each of the terms reported. The expected format for a DAVID
  result file is tab-seperated format. The following fields should be present:

       Field            -- Type   -- Example
  -----------------------------------------------------------------------------
  (0)  Category         -- string -- GOTERM_BP_FAT
  (1)  Term             -- string -- GO:0046907~intracellular transport
  (2)  Count            -- int    -- 43
  (3)  Percent          -- float  -- 11.345646437994723
  (4)  PValue           -- float  -- 1.3232857694449546E-9
  (5)  Genes            -- string -- ARSB, KPNA6, GNAS,
  (6)  List Total       -- int    -- 310
  (7)  Pop Hits	        -- int    -- 657
  (8)  Pop Total	      -- int    -- 13528
  (9)  Fold Enrichment	-- float  -- 2.8561103746256196
  (10) Bonferroni	      -- float  -- 2.6293654579179204E-6
  (11) Benjamini	      -- float  -- 2.6293654579179204E-6
  (12) FDR              -- float  -- 2.2734203852792234E-6

  The first line is a header giving the field names -- this is ignored though,
  and we expect them in the order given above.

  Most of the fields are ignored at present; we take fields 1,2, and 13 (as the
  signficiant/p-value). When parsing the term field, we try to extract a term
  ID by splitting on tilde.

  :param fn: the file to parse
  :param verbose: if True, output progress to stderr.
  """
  first = True
  for line in open(fn):
    line = line.strip()
    if line == "":
      continue
    if first:
      first = False
      continue

    parts = line.split("\t")
    if len(parts) != NUM_FIELDS_IN_DAVID_RECORD:
      raise IOError("failed to parse " + fn + " as DAVID result file. "
                    + "Expected " + str(NUM_FIELDS_IN_DAVID_RECORD) + " "
                    + "tab-separated fields, but found "
                    + str(len(parts)) + " instead")
    n_parts = parts[1].split("~")
    name = n_parts[-1].strip()
    identifier = n_parts[0] if len(n_parts) > 1 else None
    catagory = parts[0].strip()
    try:
      p_val = float(parts[PVAL_FIELD_NUM])
    except ValueError:
      raise IOError("Failed to parse " + fn + " as DAVID result file. "
                    + "Expected field " + str(PVAL_FIELD_NUM) + " "
                    + "to contain a floating point number "
                    + "(Benjamini), found this instead: "
                    + str(parts[PVAL_FIELD_NUM]))
    yield GeneOntologyEnrichmentResult(name, p_val, identifier, catagory)


###############################################################################
#                           BULK LOADING FUNCTIONS                            #
###############################################################################

def david_results_load_file(fn, verbose=False):
  """
  Load a set of DAVID gene ontology results as a list of GeneOntologyTerm
  objects

  :param fn:
  :param verbose:
  """
  return [x for x in david_results_iterator(fn, verbose)]
