================
Retrotransposons
================

Pyokit has two data structures for representing retro-transposons:
Retrotransposon and RetrotransposonOccurrence. The first of these contains
information about retro-transposons that is general to all occurrences, such
as the name, family name and concensus sequence. The second describes one
particular occurrence of a retrotransposon, including details about where
it occurs and the name of the sequence it occurrs in (usually a chromosome).
An occurrence objsect also stores the start and end co-ordinates into the
concensus sequence (becuase sometimes it is truncated at the 3' or 5' end),
and whether the match is to the concensus or the reverse compliment of the
concensus. RetrotransposonOccurrence objects contain references to
Retrotransposon objects in a one-to-many fashion. Optionally, the full
alignment between the consensus and the occurrence might also be stored
by the object, which allows more accurate liftover between the two.


-------------------
Loading annotations
-------------------
At the moment, the only annotations that are supported are those of
RepeatMasker. The coordinate annoation (which describes the locations of the
repeat occurrences), as well the full alignments can be loaded seperately,
or at the same time.

``````````````````````````````````
Repeatmasker coordinate annotation
``````````````````````````````````
Repeatmasker coordinate annotations are loaded using the following function
from the IO module. alignment_index should be left as None (default) when
loading only annotations (more details about this below).

.. autofunction:: pyokit.io.repeatMaskerIterators.repeat_masker_iterator

The description of what each field means is given below

.. autofunction:: pyokit.datastruct.retrotransposon.from_repeat_masker_string

````````````````````````````
Repeatmasker full alignments
````````````````````````````
There is an iterator in the io.alignmentIterators module for iterating over
repeatmasker alignments; you can read about how it works in
:ref:`readWriteAlignmentsSection`.

This should be sufficient for any sequential access you need to do. If you
need random access to the alignments you could use this to load them into into
whatever data-structure you like, but the large size of these files generally
makes it impractical to load the full thing. Instead, it's easier to build an
index of the file and then load the alignments you need on-demand. Building
the index is fairly straighforward and simply requires an index friendly
iterator (one which doesn't buffer the stream, and supports seeking; see the
above link on reading and writing alignments), and a hash function to
hash each element the iterator yields to a unique ID (easily done, since the
repeat-masker format contains a unique ID for each alignment)::

  from pyokit.io.indexedFile import IndexedFile
  from pyokit.io.alignmentIterators import repeat_masker_alignment_iterator

  def extract_UID(rm_alignment):
    return rm_alignment.meta[multipleAlignment.RM_ID_KEY]

  index = IndexedFile(filename, repeat_masker_alignment_iterator, extract_UID)

You can then use index to extract any alignment you want in O(1), as long as
you know the ID for it, using the subscript operator::

  index[5158]

You can read more about how indexing works in pyokit by taking a look at
:ref:`fileIndexesSection`.

This may not be all that useful on its own, since you need to know the IDs to
get the alignment you want, but it is very helpful when paired with the
repeatmasker coordinate annotations, as described in the next section.

```````````````````````````````````````````````````````
Repeatmasker annotations with on-demand full alignments
```````````````````````````````````````````````````````
The iterator for repeatmasker coordinate annotations accepts an alignment
index as a parameter. If this is provided, each RetrotransposonOccurrence
yielded by the iterator will have access to its full alignment by using the
index. This is loaded on-demand when it is needed using the index without
the application programmer needing to do anything.

--------
Liftover
--------
One operation that I often want to perform is to lift a genomic region that
overlaps a retrotransposon occurrence to the consensus sequence -- to convert
the coordinates of the region so they're relative to the consensus sequence
start. Although this is conceptually straighfoward, there are a couple of
details that complicate matters.

The first wrinkle is insertions and deletions. A retrotransposon occurrence
might have been inserted a long time ago, in which case other insertions and
deletions can occur inside of it. So although you know the start and end
coordinates of the occurrence, these may define a region that is longer or
shorter than the consensus because of insertions and deletions. Without
knowing exactly where these indels are located, your lift over will be
inaccurate. If you have the full alignment, then these can be corrected for,
but that leads to the second wrinkle: the full alignments are large; loading
them whole is generally not practical.

````````````````````````
Coordinate-only liftover
````````````````````````
The folloing program demonstrates how to compute the frequency with which
BED regions in a given file overlap each region of the consensus sequence of
retrotransposons if only the annotations are avaialable (i.e. no full
alignments)::

  import os, sys
  from pyokit.datastruct.genomicInterval import intervalTreesFromList
  from pyokit.io.bedIterators import BEDIterator
  from pyokit.io.repeatMaskerIterators import repeat_masker_iterator

  def build_profiles(repeat_masker_fn, peaks_fn):
    l = [e for e in repeat_masker_iterator(repeat_masker_fn, header=False)]
    transposon_trees = intervalTreesFromList(l, verbose=True)
    for peak in BEDIterator(peaks_fn):
      if not peak.chrom in transposon_trees : continue
      hits = transposon_trees[peak.chrom].intersectingInterval(peak.start,
                                                               peak.end)
      for hit in hits :
        sys.stderr.write("found hit for peak " + str(peak) +\
                         " in transposon " + str(hit) + "\n")
        lifted_peak = hit.liftover(peak)
        sys.stderr.write("\tconverted peak to " + str(lifted_peak) + "\n")
        for i in range(lifted_peak.start, lifted_peak.end):
          print hit.name + "\t" + hit.family_name + "\t" + str(i)

  repeat_masker_fn = sys.argv[1]
  peaks_fn = sys.argv[2]
  build_profiles(repeat_masker_fn, peaks_fn)


```````````````````````
Full-alignment liftover
```````````````````````
