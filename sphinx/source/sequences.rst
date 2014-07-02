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
   >two
   CTAGCTAGCGCGCTAGTGC
   >three
   TTTCGAGCCGGGGCAAAAA

This script will produce this output:

.. code-block:: bash

   >one
   TACGCGCTAGCGCATCAGT
   >two
   GCACTAGCGCGCTAGCTAG
   >three
   TTTTTGCCCCGGCTCGAAA

Here's the signature for the iterator:

.. autofunction:: pyokit.io.fastaIterators.fastaIterator

----------------------
Processing FastQ files
----------------------

.. autofunction:: pyokit.io.fastqIterators.fastqIterator

----------------------
Manipulating sequences
----------------------
