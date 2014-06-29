====================
Processing BED files
====================
Pyokit has a number of classes and functions for processing BED files, which are
plain-text files containing genomic intervals. There are basically two approches
for processing them in Pyokit: iterate over the BED file, which provides access
to one (or a small subset of) genomic interval at a time without storing the
whole file in memory, or bulk-load them into some data structure. Both
approaches are covered below. Pyokit also contains facilities for processing
multiple BED files simultaneously with the intention of finding matching or
missing entries in one or more of the files. Examples of this are also given
below.


You can find the definition of the BED format here:
http://genome.ucsc.edu/FAQ/FAQformat.html#format1


------------------------
Iterating over BED files
------------------------
If you only need sequential access to the genomic intervals in a BED file, then
this is the best approach -- memeory usage will be O(1).

The basic BED iterator
``````````````````````
This is the basic iterator, and will be sufficient in most situations. The
usual idiom would be:

>>> from pyokit.io.bedIterators import BEDIterator
>>> for e in BEDIterator(myFile) :
>>>   # e is a GenomicInterval object.
>>>   # Do whatever you want with it here.

Here's the signature of the basic BED iterator function:

.. autofunction:: pyokit.io.bedIterators.BEDIterator

The paired BED iterator
```````````````````````
Sometimes you may want to iterate over multiple BED files at the same time.
That's what the paired BED iterator is for. You can either ignore intervals
that don't exist in all of the BED files (the default), or you can populate
missing intervals on-the-fly with a default. By deault, intervals are considered
to match if they have the same chromosome, start and end index, and strand.
Which fields are used for matching can be adjusted though. Here's an example of
using the iterator to find all the matching intervals in a set of BED files.

>>> from itertools import izip
>>> from pyokit.io.bedIterators import pairedBEDIterator
>>>
>>> filenames = ["one.bed","two.bed","three.bed"]
>>> for lst in pairedBEDIterator(filenames) :
>>>   # each element yielded by the iterator will be a list of GenomicInterval
>>>   # objects; the index in the list will match the index of the filename that
>>>   # the interval came from. The chrom, start and end will always match
>>>   chrom = set([x.chrom for x in lst])
>>>   start = set([x.start for x in lst])
>>>   end = set([x.end for x in lst])
>>>   assert(len(chrom) == 1)
>>>   assert(len(start) == 1)
>>>   assert(len(end) == 1)
>>>   chrom = list(chrom)[0]
>>>   start = list(start)[0]
>>>   end = list(end)[0]
>>>
>>>   msgs = ", ".join([" in ".join(list(x)) for x in
>>>                     izip([x.name for x in lst], filenames)])
>>>   print "region at " + str(chrom) + " " + str(start) + " " +\
>>>         str(end) + " has name " + msgs

Given the following contents of ``one.bed``, ``two.bed`` and ``three.bed``:

.. code-block:: bash

  $ cat one.bed
  chr1	10	20	reg1	1	+
  chr1	40	50	reg2	2	-
  $ cat two.bed
  chr1	10	20	regA	1	+
  chr1	30	35	regB	5	-
  chr1	40	50	regC	2	-
  chr2	10	20	regD	2	+
  $ cat three.bed
  chr1	10	20	regW	1	+
  chr1	30	35	regX	5	-
  chr1	40	50	regY	2	-
  chr2	10	20	regZ	2	+

The output of the above script would be:

.. code-block:: bash

  region at chr1 10 20 has name reg1 in one.bed, regA in two.bed, regW in three.bed
  region at chr1 40 50 has name reg2 in one.bed, regC in two.bed, regY in three.bed

Here's the signature of the paired BED iterator function:

.. autofunction:: pyokit.io.bedIterators.pairedBEDIterator


----------------------
Bulk-loading BED files
----------------------
Sometimes you'll need to load all of the intervals in a BED file, generally
when you need random-access to them (rather than sequential). It's easy to load
them into a list (see below), but this isn't practical if you're trying to
look up intervals that fall within a certain part of the genome. An easy way to
do that is to load them into an interval tree, which then allows O(log(n)) time
to find intervals (instead of O(n) for checking a list).

.. note:: this isn't neccessarily the fastest way to do this; depending on what
          you're doing, it may be faster to sort your BED files beforehand and
          then iterate over them simultaneously. Still, interval trees are easy
          and I find they're usually fast enough in practice.

Loading into a list
```````````````````
The easiest way to do this is to use any of the above iterators and populate a
list, for example:

>>> allMyRegions = [x for x in BEDIterator(fn)]

Loading into an interval tree
`````````````````````````````
the `intervalTrees` function will load a BED file and return a dictionary of
interval trees. The keys to the dictionary are chromosome names, and each
entry is an interval tree for that chromosome. For details about the interval
tree objects, see :ref:`intervalTreesSection`.

.. autofunction:: pyokit.io.bedIterators.intervalTrees
