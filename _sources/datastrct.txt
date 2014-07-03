======================
Pyokit data-structures
======================

-----------------
Genomic Intervals
-----------------

Genomic intervals in pyokit are represented using the `GenomicInterval` class.
The basic information needed to define an interval is the chromosome it's on,
and the start and end indices for the itnerval. We consider the end interval to
be exclusive (i.e. it's not included in the interval). All intervals also have a
DNA strand, which defaults to the positive strand if not set. Additionally,
intervals can be given names and scores. The `GenomicInterval` class is
described below. For more information about using objects created from this
class, see  :ref:`genomicIntervalsSection`.

`````````````````````````
The GenomicInterval class
`````````````````````````
.. autoclass:: pyokit.datastruct.genomicInterval.GenomicInterval
   :members:


.. _intervalTreesSection:

--------------
Interval Trees
--------------
Interval trees are binary trees that allow random access lookup of intervals
that are intersected by a given point or interval.

``````````````````````
The IntervalTree class
``````````````````````

.. autoclass:: pyokit.datastruct.intervalTree.IntervalTree
   :members:

---------
Sequences
---------
Sequences in Pyokit are represented using the Sequence class, which wraps a
name and the actual sequence data, and provides a lot of the basic functionality
for manipulating sequences. Specializations exist for particular sequence
formats. At present, fasta and fastq are supported. The Sequence, FastaSequence
and FastqSequence classes are described below. For more information about
manipulating sequence data, see :ref:`sequencesSection`



.. _sequenceClassSection:

``````````````````
The Sequence class
``````````````````

.. autoclass:: pyokit.datastruct.sequence.Sequence
   :members:



.. _fastaSequenceClassSection:

```````````````````````
The FastaSequence class
```````````````````````

.. autoclass:: pyokit.datastruct.sequence.FastaSequence
   :members:



.. _fastqSequenceClassSection:

```````````````````````
The FastqSequence class
```````````````````````

.. autoclass:: pyokit.datastruct.sequence.FastqSequence
   :members:
