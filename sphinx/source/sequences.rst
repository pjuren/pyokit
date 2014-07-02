.. _sequencesSection:

====================
Biological Sequences
====================

Sequences in Pyokit are represented using the Sequence class, which wraps a
name and the actual sequence data, and provides a lot of the basic functionality
for manipulating sequences. Further specializations provide extra functionality
for fasta-format and fastq-format reads. For the class descriptions, see
:ref:`sequenceClassSection`, :ref:`fastaSequenceClassSection`, and
:ref:`fastqSequenceClassSection`

----------------------
Processing Fasta files
----------------------

Pyokit contains an iterator for fasta files. Here's an example of the usual
idiom for processing a fasta file in a scenario that requires only one pass
(in this case, reverse complement all the sequences in the fasta file).

>>> from pyokit.io.fastaIterators import fastaIterator
>>> for s in fastaIterator("sequences.fa") :
>>>   s.reverseComplement()
>>>   print s

Given the following contents for sequences.fa:

.. code-block:: bash

   >one
   ACTGATGCGCTAGCGCGTA
   CTGACGCG
   >two
   CTAGCTAGCGCGCTAGTGCGGG
   >three
   TTTCGAGCCG
   GGGCAAAAA

This script will produce this output:

.. code-block:: bash

   >one
   CGCGTCAGTACGCGCTAGC
   GCATCAGT
   >two
   CCCGCACTAGCGCGCTAGCTAG
   >three
   TTTTTGCCCC
   GGCTCGAAA

Notice that the iterator can handle sequence data split over multiple lines,
and that output formatting respects line width from the input file (with the
minor exception that line-width for a single sequence cannot be ragged, so
it takes the length of the first line for that sequence).

Here's the signature for the iterator:

.. autofunction:: pyokit.io.fastaIterators.fastaIterator

----------------------
Processing FastQ files
----------------------

This is basically the same as fasta files, but there are a few extra wrinkles.
Pyokit provides three iterators for fastq files. The difference between them
is how they handle data with sequence and/or quality data split over multiple
lines. The first, fastqIteratorSimple, ignores this possibility and expects
all fastq records to occupy 4 lines. In general this will be fine, as fastq
data, in my experience, always follows this convention. The second iterator,
fastqIteratorComplex, allows for the processing of files/streams where this
convention isn't followed (but it is slower, so should be avoided if possible).
Finally, there is a general function, called fastqIterator which at present just
wraps the simple iterator. In almost all cases, it will be sufficient just to
use this.

Here's the same example from above, but this time using a fastq iterator:

>>> from pyokit.io.fastqIterators import fastqIterator
>>> for s in fastqIterator("sequences.fq") :
>>>   s.reverseComplement()
>>>   print s

Given this input fastq file:

.. code-block:: bash

   @one
   ACTGATGCGCTAGCGCGTACTGACGCG
   +one
   !''*((((***+))%%%++)(%%%%).
   @two
   CTAGCTAGCGCGCTAGTGCGGG
   +two
   1***-+*''))**55CCF>>!*
   @three
   TTTCGAGCCGGGGCAAAAA
   +three
   CCCCCCC65+))%%%+%%)

The output will be:

.. code-block:: bash

   @one
   CGCGTCAGTACGCGCTAGCGCATCAGT
   +one
   .)%%%%()++%%%))+***((((*''!
   @two
   CCCGCACTAGCGCGCTAGCTAG
   +two
   *!>>FCC55**))''*+-***1
   @three
   TTTTTGCCCCGGCTCGAAA
   +three
   )%%+%%%))+56CCCCCCC

Notice that quality scores are also correctly reversed.

Here's the signature for the general fastQ iterator:

.. autofunction:: pyokit.io.fastqIterators.fastqIterator

Here's the signature for the simple iterator:

.. autofunction:: pyokit.io.fastqIterators.fastqIteratorSimple

And here's the signature for the iterator that handles multi-line quality and
sequence data:

.. autofunction:: pyokit.io.fastqIterators.fastqIteratorComplex

----------------------
Manipulating sequences
----------------------

Pretty much everything you need should be covered in the class descriptions for
the sequence classes: see :ref:`sequenceClassSection`,
:ref:`fastaSequenceClassSection`, and :ref:`fastqSequenceClassSection`
