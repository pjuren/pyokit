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


The basic BED iterator
``````````````````````
This is the basic iterator, and will be sufficient in most situations. The
usual idiom would be

>> from pyokit.mapping.bedIterators import BEDIterator
>> for e in BEDIterator(myFile) :
>>   # e is a GenomicInterval object.
>>   # Do whatever you want with it here.

.. autofunction:: pyokit.mapping.bedIterators.BEDIterator



----------------------
Bulk-loading BED files
----------------------
Loading into a list
```````````````````
The easiest way to do this is to use any of the above iterators and populate a
list, for example:

Loading into an interval tree
`````````````````````````````
.. autofunction:: pyokit.mapping.bedIterators.intervalTrees
