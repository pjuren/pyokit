==============================
Manipulating Genomic Intervals
==============================

Pyokit has a class for representing and manipulating genomic intervals, and some
associated functions for manipulating collections of genomic intervals.

------------------------
Collapsing interval sets
------------------------
.. autofunction:: pyokit.datastruct.genomicInterval.collapseRegions

--------------------------
Intersecting interval sets
--------------------------
.. autofunction:: pyokit.datastruct.genomicInterval.regionsIntersection


-----------------------
Searching interval sets
-----------------------
Usually when you're searching for particular intervals, you're looking for
ones that exist in some region of the genome, or contain some point in the
genome. The naive approach would be to check all intervals, requiring O(n) time.
Although this sounds inefficient, it might be the best option if you only
need to make one pass of the intervals. Here's an example:

>>> from pyokit.io.bedIterators import BEDIterator
>>>
>>> # print out all intervals that intersect chr1, 10, 20
>>> roi = GenomicInterval("chr1", 10, 20)
>>> hits = [e for e in BEDIterator("somefile.bed") if e.intersects(roi)]

More often though, you'll want to make multiple passes. For example, you might
want to count how many intervals in ``two.bed`` intersect each interval in
``one.bed``. The naive approach is to use nested iterators, but it has
worst case complexity of O(n^2). This occurs sufficiently often that Pyokit has
an iterator for dealing with exactly this situation: the bucket iterator --
which will iterate over one set of genomic intervals (the buckets), and for each
one will give you the list of intervals in another set that are intersected by
it. You can provide it with regular BEDIterators, or lists, or any iterable
object that contains genomic intervals. Here's the signature of the iterator:

.. autofunction:: pyokit.datastruct.genomicInterval.bucketIterator

As an example of using this iterator, consider that you have the following
sets of buckets, and elements that you want to bin:

>>> # set up some buckets
>>> g1 = GenomicInterval("chr1", 1,  6)
>>> g2 = GenomicInterval("chr1", 4,  10)
>>> g3 = GenomicInterval("chr1", 5,  19)
>>> g4 = GenomicInterval("chr1", 9,  15)
>>> g5 = GenomicInterval("chr1", 18, 22)
>>> g6 = GenomicInterval("chr2", 1,  12)
>>> g7 = GenomicInterval("chr2", 4,  7)
>>> buckets = [g1,g2,g3,g4,g5,g6,g7]
>>>
>>> # set up some elements to bin
>>> e1 = GenomicInterval("chr1", 2,  3)
>>> e2 = GenomicInterval("chr1", 4,  5)
>>> e3 = GenomicInterval("chr1", 4,  5)
>>> e4 = GenomicInterval("chr1", 7,  9)
>>> e5 = GenomicInterval("chr1", 7,  8)
>>> e6 = GenomicInterval("chr1", 10, 14)
>>> e7 = GenomicInterval("chr1", 12, 20)
>>> e8 = GenomicInterval("chr1", 20, 23)
>>> e9 = GenomicInterval("chr1", 24, 26)
>>> e10 = GenomicInterval("chr2", 1,  3)
>>> e11 = GenomicInterval("chr2", 2,  6)
>>> e12 = GenomicInterval("chr2", 4,  8)
>>> e13 = GenomicInterval("chr2", 10, 15)
>>> elements = [e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13]

Now we can bin these elements as follows:

>>> actual = [(x, set(l)) for x,l in bucketIterator(elements, buckets)]
>>> actualStr = "\n".join([str(x) + " == " + ", ".join([str(y) for y in l])
>>>                        for x, l in actual])
>>> print actualStr

This will give you the following:

.. code-block:: bash

   chr1	1	6 == chr1	4	5, chr1	2	3
   chr1	4	10 == chr1	4	5, chr1	7	9,chr1	7	8
   chr1	5	19 == chr1	10	14, chr1	12	20, chr1	7	9, chr1	7	8
   chr1	9	15 == chr1	12	20, chr1	10	14
   chr1	18	22 == chr1	20	23, chr1	12	20
   chr2	1	12 == chr2	2	6, chr2	4	8, chr2	10	15, chr2	1	3
   chr2	4	7 == chr2	2	6, chr2	4	8


In certain more complicated scenarios, it might be easier to use an interval
tree to do this kind of thing. Here is code using the interval tree from list
function to perform the same operation:

>>> trees = intervalTreesFromList(elements, openEnded = True)
>>> treeRes = [(g, set(trees[g.chrom].intersectingInterval(g.start, g.end)))
>>>            for g in buckets]
>>> treeReStr = "\n".join([str(x) + " == " + ",".join([str(y) for y in l])
>>>                        for x, l in treeRes])
>>> print treeReStr

The output will be the same.
