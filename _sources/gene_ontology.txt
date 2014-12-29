=============
Gene Ontology
=============
Gene ontology data is represented by two data types; you can find details
of these in :ref:`ontologyClassSection`. At present, the only input
format supported is that produced by DAVID. Here's an example of loading
several DAVID output files and producing a new data frame containing the
GO terms, catagories, and significance (p-value) from each enrichment analysis
for those terms where at least one of the analyses was significant:

>>> import os, sys
>>> from pyokit.io.david import david_results_iterator
>>> PVAL_THRESHOLD = 0.01
>>> filenames = sys.argv[1:]
>>> # load all of the DAVID results for each file
>>> by_trm = {}
>>> for fn in filenames:
>>>   for r in david_results_iterator(fn):
>>>     if not r.name in by_trm:
>>>       by_trm[r.name] = {}
>>>     by_trm[r.name][fn] = r
>>> # drop terms where no file has p < threshold
>>> by_trm = {term:by_trm[term] for term in by_trm
>>>           if min([by_term[term][fn].pvalue for fn in by_term[term]]) < PVAL_THRESHOLD}
>>> # output
>>> for term in by_trm:
>>>   for fn in by_trm[term]:
>>>     r = by_trm[term][fn]
>>>     print r.name + "\t" + str(r.pvalue) + "\t" + r.catagory + "\t" + fn

This makes use of an iterator for DAVID otuput-format files. Here are the
details of that function:

.. autofunction:: pyokit.io.david.david_results_iterator
