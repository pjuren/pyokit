#!/usr/bin/python

"""
Date of Creation: 3rd Apr 2010.

Description:   Classes and functions for manipulating Genomic Intervals

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

# standard python library imports
import sys
import unittest
import copy
from heapq import heappush, heappop

# pyokit imports
from pyokit.util.progressIndicator import ProgressIndicator
from pyokit.datastruct.intervalTree import IntervalTree
from pyokit.util.meta import AutoApplyIterator


###############################################################################
#                            EXCEPTION CLASSES                                #
###############################################################################

class GenomicIntervalError(Exception):

  """Errors raised when manipulating genomic intervals."""

  def __init__(self, msg):
    """Constructor for genomic interval errors."""
    self.value = msg

  def __str__(self):
    """Get string representation of this exception."""
    return repr(self.value)


###############################################################################
#        FUNCTIONS FOR MANIPULATING COLLECTIONS OF GENOMIC INTERVALS          #
###############################################################################

def intervalTreesFromList(inElements, verbose=False, openEnded=False):
  """
  build a dictionary, indexed by chrom name, of interval trees for each chrom.

  :param inElements: list of genomic intervals. Members of the list must have
                     chrom, start and end fields; no other restrictions.
  :param verbose: output progress messages to sys.stderr if True
  """
  elements = {}
  if verbose:
    totalLines = len(inElements)
    pind = ProgressIndicator(totalToDo=totalLines,
                             messagePrefix="completed",
                             messageSuffix="of parsing")

  for element in inElements:
    if element.chrom not in elements:
      elements[element.chrom] = []
    elements[element.chrom].append(element)
    if verbose:
      pind.done += 1
      pind.showProgress()

  # create an interval tree for each list
  trees = {}
  if verbose:
    totalLines = len(elements)
    pind = ProgressIndicator(totalToDo=totalLines,
                             messagePrefix="completed",
                             messageSuffix="of making interval trees")
  for chrom in elements:
    trees[chrom] = IntervalTree(elements[chrom], openEnded)
    if verbose:
      pind.done += 1
      pind.showProgress()

  return trees


def collapseRegions(s):
  """
  Get the union of a set of genomic intervals.

  given a list of genomic intervals with chromosome, start and end field,
  collapse into a set of non-overlapping intervals. Intervals must be sorted by
  chromosome and then start coordinate.

  :note: O(n) time, O(n) space
  :return:  list of intervals that define the collapsed regions. Note that
            these are all new objects, no existing object from s is returned
            or altered. Returned regions will all have name "X", strand +
            and score 0
  :param s: list of genomic regions to collapse
  :raise GenomicIntervalError: if the input regions are not correctly sorted
                               (chromosome then start)
  """
  debug = False

  if len(s) == 0 or len(s) == 1:
    return copy.deepcopy(s)

  res = []
  current = copy.copy(s[0])
  current.strand = '+'
  current.score = 0
  current.name = "X"
  for i in range(1, len(s)):
    if debug:
      sys.stderr.write("processing " + str(s[i]) + "\n")
      sys.stderr.write("\tcurrent is: " + str(current) + "\n")

    # make sure things are sorted..
    if (s[i].chrom < s[i - 1].chrom) or \
       (s[i].chrom == s[i - 1].chrom and s[i].start < s[i - 1].start):
      raise GenomicIntervalError("collapsing regions failed. saw this " +
                                 "region: " + str(s[i - 1]) + " before this " +
                                 "one: " + str(s[i]))

    # because of sorting order, we know that nothing else exists with
    # start less than s[i] which we haven't already seen.
    if s[i].start > current.end or s[i].chrom != current.chrom:
      res.append(current)
      current = copy.copy(s[i])
      current.strand = '+'
      current.score = 0
      current.name = "X"
    else:
      current.end = max(s[i].end, current.end)

  # don't forget the last one...
  res.append(current)

  return res


def regionsIntersection(s1, s2):
  """
  given two lists of genomic regions with chromosome, start and end
  coordinates, return a new list of regions which is the intersection of those
  two sets. Lists must be sorted by chromosome and start index

  :return: new list that represents the intersection of the two input lists.
           output regions will all have name "X", be on strand "+" and have
           score 0
  :param s1: first list of genomic regions
  :param s2: second list of genomic regions
  :raise GenomicIntervalError: if the input regions are not sorted correctly
                               (by chromosome and start index)
  :note: O(n) time, O(n) space; informally, might use up to 3x space of input
  """
  debug = False

  # we don't need to explicitly check for sorting because sorted order is
  # a post-condition of the collapsing function
  s1_c = collapseRegions(s1)
  s2_c = collapseRegions(s2)

  if len(s1_c) == 0 or len(s2_c) == 0:
    return []

  res = []
  j = 0
  for i in range(0, len(s1_c)):
    if debug:
      sys.stderr.write("processing from s1_c : " + str(s1_c[i]) + "\n")

    # find first thing in s2_c with end in or after s1_c[i]
    if debug:
      sys.stderr.write("i = " + str(i) + " and j = " + str(j) + "\n")
    while (j < len(s2_c) and
           (s2_c[j].chrom < s1_c[i].chrom or
           (s2_c[j].chrom == s1_c[i].chrom and s2_c[j].end <= s1_c[i].start))):
      j += 1
    # nothing intersects if we hit the end of s2, or the end of the chrom,
    # or we're still on the same chrom but start after the end of s2_c[i]
    if j >= len(s2_c) or s2_c[j].chrom > s1_c[i].chrom or \
       (s2_c[j].chrom == s1_c[i].chrom and s2_c[j].start >= s1_c[i].end):
      continue

    # now everything at or after j in s2_c that starts before
    # the end of s1_c must overlap with it
    while s2_c[j].start < s1_c[i].end:
      s = max(s1_c[i].start, s2_c[j].start)
      e = min(s1_c[i].end, s2_c[j].end)
      overlap = GenomicInterval(s1_c[i].chrom, s, e, "X", 0, "+")
      if debug:
        sys.stderr.write("\tadding to overlaps: " + str(overlap) + "\n")
      res.append(overlap)
      j += 1
      if j >= len(s2_c) or s2_c[j].chrom != s1_c[i].chrom:
        break

    # it's possible the last intersecting element runs on to the
    # next element from s1_c, so...
    j -= 1
    if debug:
      sys.stderr.write("\tmoving s2_c index back to " + str(s2_c[j]) + "\n")

  return res


def bucketIterator(elements, buckets):
  """
  For each bucket in buckets, yield it and any elements that overlap it.

  :param elements: the genomic intervals to place into the buckets. Must be
                   sorted by chromosome and start index. This could be a
                   list, or an iterator.
  :param buckets: the buckets into which genomic intervals should be binned.
                  Must be sorted by chromosome and start index. This could be
                  a list, or an iterator
  :return: iterator that will yeild a tuple of 1 bucket and 1 list of elements
           in the bucket for each call to __next__().
  """
  def check_sorted(current, previous):
    if (previous is not None) and \
       ((previous.chrom > current.chrom) or
        ((previous.chrom == current.chrom) and
         (previous.start > current.start))):
      raise GenomicIntervalError("not sorted")

  def updateOpen(openHeap, elementIterator, bucketChrom,
                 bucketStart, bucketEnd):
    """
    Drop elements from heap which start earlier than current bucket.

    Update the open heap so that it contains only elements that end after the
    start of the current bucket. Note that the heap may already contain some
    elements that start after the end of the current bucket, if a previous
    bucket ended after the end of this one and brought them into the set.

    :param openHeap: a min heap of elements; uses the default sorting order
                     for the genomic intervals, which is by end index. This
                     is what we're updating.
    :param elementIterator: an iterator from which we will pull new elements.
                            Elements yielded by this iterator must be sorted
                            by start index. Must be 'peakable'
    :param bucketChrom: the chromosome of the current bucket.
    :param bucketStart: the start index of the current bucket.
    :param bucketEnd: the end index of the current bucket.
    """
    # first, we're going to pop elements from the heap which can no longer
    # overalp this or any future buckets. Buckets are sorted by start, so
    # we'll never see another bucket that starts earlier than this one --
    # hence any elements that end before the start of this bucket will never be
    # used again and can be dropped. Elements in the heap are ordered by end
    # index, so once we reach an element in the heap that does not end before
    # the start of this bucket, we are sure that no others will come after it
    # which do end before the start of this bucket. So we can stop dropping.
    while len(openHeap) > 0 and ((openHeap[0].chrom < bucketChrom) or
                                 ((openHeap[0].chrom == bucketChrom) and
                                 (openHeap[0].end <= bucketStart))):
      heappop(openHeap)

    # now we're going to add new elements from the iterator to the heap. As
    # we know that elements in the iterator are sorted by start index, we know
    # that once we see an element that has a start index greater than the end
    # of this bucket, we can stop -- everything else after it will also start
    # after the end of this bucket.
    while (elementIterator.peek() != None) and \
          ((elementIterator.peek().chrom < bucketChrom) or
           ((elementIterator.peek().chrom == bucketChrom) and
            (elementIterator.peek().start < bucketEnd))):
      e = elementIterator.__next__()
      # if e falls before this bucket, we can skip it; buckets are sorted by
      # start, so no other buckets start earlier than this one and so it
      # cannot intersect any others.
      if (e.chrom < bucketChrom) or \
         (e.chrom == bucketChrom and e.end <= bucketStart):
         continue
      # now we know e intersects this bucket..
      heappush(openHeap, e)

  openElems = []
  prevBucket = None
  elementIter = AutoApplyIterator(elements, check_sorted)
  for bucket in buckets:
    # make sure the buckets are sorted by start index
    if prevBucket is not None and ((bucket.chrom < prevBucket.chrom) or
                                   (bucket.chrom == prevBucket.chrom and
                                    bucket.start < prevBucket.start)):
      raise GenomicIntervalError("not sorted")
    updateOpen(openElems, elementIter, bucket.chrom, bucket. start, bucket.end)

    # be careful here not to leak a reference to the heap; if the caller
    # decides to mess with that list, it'll screw us up. Anyway, we need a
    # final check here to make sure we trim off any elements that exceed the
    # end of this bucket.
    yield bucket, [x for x in openElems if x.start < bucket.end]
    prevBucket = bucket


###############################################################################
#                 FUNCTIONS FOR PARSING GENOMIC INTERVALS                     #
###############################################################################

def parseWigString(line, scoreType=int):
  """
  Parse a string in simple Wig format and return a GenomicInterval.

  :param line: the string to be parsed
  :param scoreType: treat the score field as having this type.
  :return: GenomicInterval object representing this wig line; the name of the
           interval will be set to 'X', and it's strand to the default.
  """
  parts = line.split("\t")
  if (len(parts) < 4):
    raise GenomicIntervalError("failed to parse " + line +
                               " as wig format, too few fields")
  return GenomicInterval(parts[0].strip(), int(parts[1]), int(parts[2]), None,
                         scoreType(parts[3]))


def parseBEDString(line, scoreType=int, dropAfter=None):
  """
  Parse a string in BED format and return a GenomicInterval object.

  :param line: the string to be parsed
  :param dropAfter: an int indicating that any fields after and including
                    this field should be ignored as they don't conform to
                    the BED format. By default, None, meaning we use all
                    fields. Index from zero.
  :return: GenomicInterval object built from the BED string representation
  """
  peices = line.split("\t")
  if dropAfter is not None:
    peices = peices[0:dropAfter]
  if len(peices) < 3:
    raise GenomicIntervalError("BED elements must have at least chrom, " +
                               "start and end; found only " +
                               str(len(peices)) + " in " + line)
  chrom = peices[0]
  start = peices[1]
  end = peices[2]

  name = None
  score = None
  strand = None

  if len(peices) >= 4 is not None:
    name = peices[3]
  if len(peices) >= 5 is not None:
    score = peices[4]
  if len(peices) >= 6 is not None:
    strand = peices[5]

  return GenomicInterval(chrom, start, end, name, score, strand, scoreType)


###############################################################################
#                           GENOMIC INTERVAL CLASS                            #
###############################################################################

class GenomicInterval(object):

  """
  Represents contiguous segment of a genome; inclusive of start, but not end.

  :param chrom: The chromosome this genomic interval is on.
  :param start: The start of this genomic interval (inclusive).
  :param end: The end of this genomic interval (exclusive).
  :param name: A name associated with this genomic interval.
  :param score: A score associated with this genomic interval.
  :param strand: The DNA strand (+ or -) that this genomic interval is on.
  :param scoreType: The type (e.g. int, float) of the score associated with
                    this genomic interval
  """

  POSITIVE_STRAND = "+"
  NEGATIVE_STRAND = "-"
  DEFAULT_STRAND = POSITIVE_STRAND

  def __init__(self, chrom, start, end, name=None, score=None,
               strand=None, scoreType=int):
    """Constructor for the GenomicInterval. See class docstring for details."""
    # the basic info -- we must get at least this much
    if chrom is None or start is None or end is None:
      raise GenomicIntervalError("Must provided at least chrom, start, end " +
                                 "for Genomic Interval")
    self.chrom = chrom.strip()
    self.start = int(start)
    self.end = int(end)

    # we might get the following too...
    # name:
    self.name = None
    if name is not None:
      self.name = name.strip()

    # score
    self.score = None
    if score is not None:
      self.score = scoreType(score)

    # strand
    self.strand = None
    if strand is not None:
      self.strand = strand.strip()
    if self.strand is not None and self.strand != self.POSITIVE_STRAND and \
       self.strand != self.NEGATIVE_STRAND:
      raise GenomicIntervalError("Invalid strand: " + self.strand +
                                 "; expected either " + self.POSITIVE_STRAND +
                                 " or " + self.NEGATIVE_STRAND)

  def __hash__(self):
    """:return: a hash of this GenomicInterval."""
    return hash(str(self))

  def __eq__(self, e):
    """:return: True if self is equal to e."""
    if e is None:
      return False
    try:
      return (self.chrom == e.chrom and self.start == e.start and
              self.end == e.end and self.name == e.name and
              self.score == e.score and self.strand == e.strand)
    except AttributeError:
      return False

  def __lt__(self, rhs):
    """
    Check whether self < rhs? Comparison is first by chrom, then by end index.

    :param rhs: object from the right-hand-side of the comparaison self < rhs
    :return: True if self is considered 'less' than rhs.
    """
    if rhs is None:
      return False
    if self.chrom < rhs.chrom:
      return True
    if self.chrom > rhs.chrom:
      return False
    if self.end < rhs.end:
      return True
    return False

  def sameRegion(self, e):
    """
    Check whether self represents the same DNA region as e.

    :param e: genomic region to compare against
    :return: True if self and e are for the same region (ignores differences
             in non-region related fields, such as name or score -- but does
             consider strand)
    """
    if e is None:
      return False
    return (self.chrom == e.chrom and self.start == e.start and
            self.end == e.end and self.name == e.name and
            self.strand == e.strand)

  def __len__(self):
    """:return the length of this interval; number of nucleotides in it."""
    return (self.end - self.start)

  def __str__(self):
    """
    Produce a string representation of the interval.

    Only those fields which we have values for will be output.

    :return: String representation of the genomic interval in BED format.
    """
    delim = "\t"
    res = self.chrom + delim + str(self.start) + delim + str(self.end)
    if self.name is not None:
      res += (delim + str(self.name))
    if self.score is not None:
      res += (delim + str(self.score))
    if self.strand is not None:
      res += (delim + str(self.strand))

    return res

  def distance(self, e):
    """
    Distance between this interval and e -- number of nucleotides.

    We consider intervals that overlap to have a distance of 0 to each other.
    The distance between two intervals on different chromosomes is considered
    undefined, and causes an exception to be raised.

    :return: the distance from this GenomicInterval to e.
    :param e: the other genomic interval to find the distance to.
    :raise GenomicIntervalError: if self and e are on different chromosomes.
    """
    if e.chrom != self.chrom:
      raise GenomicIntervalError("cannot get distance from " + str(self) +
                                 " to " + str(e) + " as they're on " +
                                 "different chromosomes")
    dist = 0
    if not e.intersects(self):
      dist = max(self.start, e.start) - min(self.end, e.end)
    return dist

  def signedDistance(self, e):
    """
    Signed distance between this interval and e -- number of nucleotides.

    Get the signed distance from this genomic interval to another one. We
    consider intervals that overlap to have a distance of 0 to each other. The
    distance between two intervals on different chromosomes is considered
    undefined, and causes an exception to be raised. If e comes earlier than
    self, the distance will be negative.

    :return: the signed distance from this GenomicInterval to e.
    :param e: the other genomic interval to find the distance to.
    :raise GenomicIntervalError: if self and e are on different chromosomes.
    """
    dist = self.distance(e)
    if e < self:
      dist = dist * -1
    return dist

  def __singleIntersect(self, e):
    """
    Get a list of GenomicRegions produced by subtracting e from self.

    :return: a list of regions that represent the result of subtracting e from
             self. The list might be empty if self is totally contained in e,
             or may contain two elements if e is totally contained in self.
             otherwise there'll be one element.
    """
    if e.chrom != self.chrom or e.end < self.start or e.start > self.end:
      # no intersection
      return [copy.copy(self)]
    if e.start <= self.start and e.end >= self.end:
      # whole self is removed.
      return []
    if e.start > self.start and e.end < self.end:
      # splits self in two
      r1 = copy.copy(self)
      r2 = copy.copy(self)
      r1.end = e.start
      r2.start = e.end
      return [r1, r2]
    if e.start <= self.start:
      # cuts off start of self
      r = copy.copy(self)
      r.start = e.end
      return [r]
    if e.end >= self.end:
      # cuts off end of self
      r = copy.copy(self)
      r.end = e.start
      return [r]
    # oops, we screwed up if we got to here..
    raise GenomicIntervalError("fatal error - failed BED subtraction of " +
                               str(e) + " from " + str(self))

  def subtract(self, es):
    """
    Subtract the BED elements in es from self.

    :param es: a list of BED elements (or anything with chrom, start, end)
    :return: a list of BED elements which represent what is left of
             self after the subtraction. This might be an empty list.
    """
    workingSet = [self]
    for e in es:
      newWorkingSet = []
      for w in workingSet:
        newWorkingSet += w.__singleIntersect(e)
      workingSet = newWorkingSet
    return workingSet

  def sizeOfOverlap(self, e):
    """
    Get the size of the overlap between self and e.

    :return: the number of bases that are shared in common between self and e.
    """
    # no overlap
    if not self.intersects(e):
      return 0

    # complete inclusion..
    if e.start >= self.start and e.end <= self.end:
      return len(e)
    if self.start >= e.start and self.end <= e.end:
      return len(self)

    # partial overlap
    if e.start > self.start:
      return (self.end - e.start)
    if self.start > e.start:
      return (e.end - self.start)

  def intersects(self, e):
    """
    Check whether e intersects self.

    :return: true if this elements intersects the element e
    """
    if self.chrom != e.chrom:
      return False

    if e.start >= self.start and e.start < self.end:
      return True
    if e.end > self.start and e.end <= self.end:
      return True
    if e.start < self.start and e.end > self.end:
      return True
    return False

  def isPositiveStrand(self):
    """
    Check if this genomic region is on the positive strand.

    :return: True if this element is on the positive strand
    """
    if self.strand is None and self.DEFAULT_STRAND == self.POSITIVE_STRAND:
      return True
    return self.strand == self.POSITIVE_STRAND

  def isNegativeStrand(self):
    """
    Check if this genomic interval is on the negative strand.

    :return: True if this element is on the negative strand
    """
    if self.strand is None and self.DEFAULT_STRAND == self.NEGATIVE_STRAND:
      return True
    return self.strand == self.NEGATIVE_STRAND

  def transform_center(self, size):
    """
    Tranform self so it is centered on the same spot, but has new size.

    If the region grows, the extra nucleotides will be distributed evenly to
    the 5' and 3' side of the current region if possible. If the extra is odd,
    the 3' side will get the extra one. Similarly, if the resize shrinks the
    interval, bases will be removed from the 5' and 3' sides equally; if the
    number to remove is odd, the extra one will be removed from the 3' side.

    :param size: size of the region after transformation.
    """
    if size < 1:
      raise GenomicIntervalError("Cannot resize genomic interval to " +
                                 str(size) + " bases; must be at least 1 " +
                                 "base in length")
    if size == len(self):
      return
    elif size > len(self):
      extra = size - len(self)
      extra_5_prime = extra / 2
      extra_3_prime = extra / 2 if extra % 2 == 0 else (extra / 2) + 1
      assert(extra_5_prime + extra_3_prime == extra)
      self.start = (self.start - extra_5_prime if self.isPositiveStrand()
                    else self.start - extra_3_prime)
      self.end = (self.end + extra_3_prime if self.isPositiveStrand()
                  else self.end + extra_5_prime)
    else:
      less = len(self) - size
      less_5_prime = less / 2
      less_3_prime = less / 2 if less % 2 == 0 else (less / 2) + 1
      assert(less_5_prime + less_3_prime == less)
      self.start = (self.start + less_5_prime if self.isPositiveStrand()
                    else self.start + less_3_prime)
      self.end = (self.end - less_3_prime if self.isPositiveStrand()
                  else self.end - less_5_prime)


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestGenomicInterval(unittest.TestCase):

  """Unit tests for functions and classes in this module."""

  def testCollapse(self):
    """Test collapsing regions to get union."""
    debug = False
    elements = ["chr1" + "\t" + "10" + "\t" + "20" + "\t" + "R01" + "\t" +
                "0" + "\t" + "+",
                "chr1" + "\t" + "30" + "\t" + "40" + "\t" + "R02" + "\t" +
                "1" + "\t" + "+",
                "chr1" + "\t" + "35" + "\t" + "50" + "\t" + "R03" + "\t" +
                "0" + "\t" + "+",
                "chr1" + "\t" + "45" + "\t" + "65" + "\t" + "R04" + "\t" +
                "0" + "\t" + "+",
                "chr1" + "\t" + "55" + "\t" + "60" + "\t" + "R05" + "\t" +
                "3" + "\t" + "-",
                "chr1" + "\t" + "70" + "\t" + "80" + "\t" + "R06" + "\t" +
                "0" + "\t" + "+",
                "chr1" + "\t" + "75" + "\t" + "95" + "\t" + "R07" + "\t" +
                "0" + "\t" + "+",
                "chr1" + "\t" + "85" + "\t" + "90" + "\t" + "R08" + "\t" +
                "1" + "\t" + "-",
                "chr2" + "\t" + "40" + "\t" + "60" + "\t" + "R10" + "\t" +
                "0" + "\t" + "+",
                "chr3" + "\t" + "10" + "\t" + "20" + "\t" + "R11" + "\t" +
                "4" + "\t" + "+",
                "chr3" + "\t" + "20" + "\t" + "30" + "\t" + "R12" + "\t" +
                "0" + "\t" + "-"]
    expect_str = ["chr1" + "\t" + "10" + "\t" + "20" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr1" + "\t" + "30" + "\t" + "65" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr1" + "\t" + "70" + "\t" + "95" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr2" + "\t" + "40" + "\t" + "60" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr3" + "\t" + "10" + "\t" + "30" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+"]
    input_x = [parseBEDString(x) for x in elements]
    expect = [parseBEDString(x) for x in expect_str]
    got = collapseRegions(input_x)
    if debug:
      sys.stderr.write("expect:\n")
      sys.stderr.write("\n".join([str(x) for x in expect]) + "\n")
      sys.stderr.write("got:\n")
      sys.stderr.write("\n".join([str(x) for x in got]) + "\n")
    self.assertEqual(expect, got)

  def testRegionsIntersection(self):
    """Test getting the intersection of a set of regions."""
    debug = False
    s1_elements = ["chr1" + "\t" + "40" + "\t" + "90" + "\t" + "R11" + "\t" +
                   "0" + "\t" + "+",
                   "chr1" + "\t" + "100" + "\t" + "120" + "\t" + "R12" +
                   "\t" + "3" + "\t" + "+",
                   "chr1" + "\t" + "160" + "\t" + "190" + "\t" + "R13" +
                   "\t" + "1" + "\t" + "-",
                   "chr1" + "\t" + "200" + "\t" + "210" + "\t" + "R14" +
                   "\t" + "0" + "\t" + "-",
                   "chr2" + "\t" + "10" + "\t" + "20" + "\t" + "R15" + "\t" +
                   "0" + "\t" + "-",
                   "chr3" + "\t" + "10" + "\t" + "80" + "\t" + "R16" + "\t" +
                   "1" + "\t" + "+",
                   "chr4" + "\t" + "20" + "\t" + "30" + "\t" + "R17" + "\t" +
                   "1" + "\t" + "+",
                   "chr4" + "\t" + "40" + "\t" + "50" + "\t" + "R18" + "\t" +
                   "1" + "\t" + "-",
                   "chr5" + "\t" + "40" + "\t" + "50" + "\t" + "R19" + "\t" +
                   "1" + "\t" + "-"]
    s2_elements = ["chr1" + "\t" + "10" + "\t" + "20" + "\t" + "R21" + "\t" +
                   "0" + "\t" + "+",
                   "chr1" + "\t" + "30" + "\t" + "50" + "\t" + "R22" + "\t" +
                   "1" + "\t" + "-",
                   "chr1" + "\t" + "60" + "\t" + "70" + "\t" + "R23" + "\t" +
                   "0" + "\t" + "+",
                   "chr1" + "\t" + "80" + "\t" + "110" + "\t" + "R24" + "\t" +
                   "4" + "\t" + "+",
                   "chr1" + "\t" + "130" + "\t" + "140" + "\t" + "R25" + "\t" +
                   "0" + "\t" + "+",
                   "chr1" + "\t" + "150" + "\t" + "170" + "\t" + "R26" + "\t" +
                   "1" + "\t" + "-",
                   "chr1" + "\t" + "180" + "\t" + "220" + "\t" + "R27" + "\t" +
                   "0" + "\t" + "-",
                   "chr2" + "\t" + "30" + "\t" + "40" + "\t" + "R28" + "\t" +
                   "0" + "\t" + "-",
                   "chr3" + "\t" + "20" + "\t" + "30" + "\t" + "R29" + "\t" +
                   "0" + "\t" + "-",
                   "chr3" + "\t" + "40" + "\t" + "50" + "\t" + "R210" + "\t" +
                   "0" + "\t" + "-",
                   "chr3" + "\t" + "60" + "\t" + "70" + "\t" + "R211" + "\t" +
                   "0" + "\t" + "+",
                   "chr4" + "\t" + "10" + "\t" + "60" + "\t" + "R212" + "\t" +
                   "0" + "\t" + "+",
                   "chr5" + "\t" + "10" + "\t" + "20" + "\t" + "R213" + "\t" +
                   "1" + "\t" + "-"]
    expect_str = ["chr1" + "\t" + "40" + "\t" + "50" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr1" + "\t" + "60" + "\t" + "70" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr1" + "\t" + "80" + "\t" + "90" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr1" + "\t" + "100" + "\t" + "110" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr1" + "\t" + "160" + "\t" + "170" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr1" + "\t" + "180" + "\t" + "190" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr1" + "\t" + "200" + "\t" + "210" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr3" + "\t" + "20" + "\t" + "30" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr3" + "\t" + "40" + "\t" + "50" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr3" + "\t" + "60" + "\t" + "70" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr4" + "\t" + "20" + "\t" + "30" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+",
                  "chr4" + "\t" + "40" + "\t" + "50" + "\t" + "X" + "\t" +
                  "0" + "\t" + "+"]
    input_s1 = [parseBEDString(x) for x in s1_elements]
    input_s2 = [parseBEDString(x) for x in s2_elements]
    expect = [parseBEDString(x) for x in expect_str]
    got = regionsIntersection(input_s1, input_s2)
    if debug:
      sys.stderr.write("expect:\n")
      sys.stderr.write("\n".join([str(x) for x in expect]) + "\n")
      sys.stderr.write("got:\n")
      sys.stderr.write("\n".join([str(x) for x in got]) + "\n")
    self.assertEqual(expect, got)

  def testSizeOfOverlap(self):
    """Test getting the sie of overlap between two regions."""
    # region 2 entirely inside region 1 -- ans = size of region 2
    region1 = GenomicInterval("chr1", 10, 20, "read", 1, "+")
    region2 = GenomicInterval("chr1", 10, 20, "read", 1, "+")
    self.assertEqual(region1.sizeOfOverlap(region2), len(region2))
    self.assertEqual(region1.sizeOfOverlap(region2),
                     region2.sizeOfOverlap(region1))

  def testSubtract(self):
    """Test subtracting one region from another."""
    debug = False
    a = GenomicInterval("chr1", 10, 100)
    b1 = GenomicInterval("chr1", 1, 8)
    b2 = GenomicInterval("chr1", 15, 20)
    b3 = GenomicInterval("chr1", 50, 60)
    b4 = GenomicInterval("chr1", 90, 120)
    c1 = GenomicInterval("chr2", 30, 40)

    res = a.subtract([b1, b2, b3, b4, c1])
    expect = [GenomicInterval("chr1", 10, 15),
              GenomicInterval("chr1", 20, 50),
              GenomicInterval("chr1", 60, 90)]
    res.sort()
    expect.sort()
    if debug:
      sys.stderr.write("\nDebugging BED subtraction test\n")
      sys.stderr.write("expect\n" + "\n".join([str(x) for x in expect]) + "\n")
      sys.stderr.write("got\n" + "\n".join([str(x) for x in res]) + "\n")
    assert(res == expect)

  def testBucketIterator(self):
    """
    Test the bucket iterator.

    chr 1 'genes':
    1   5   9    14  18 21
    -----   ------
       ------        ----
        --------------
     * *  **   ********
       *  *  ****      *** **

    chr 2 'genes'
    1         11
    -----------
       ---
    ** ****
     ****    *****
    """
    debug = False

    # set up some 'genes'
    g1 = GenomicInterval("chr1", 1, 6)
    g2 = GenomicInterval("chr1", 4, 10)
    g3 = GenomicInterval("chr1", 5, 19)
    g4 = GenomicInterval("chr1", 9, 15)
    g5 = GenomicInterval("chr1", 18, 22)
    g6 = GenomicInterval("chr2", 1, 12)
    g7 = GenomicInterval("chr2", 4, 7)
    buckets = [g1, g2, g3, g4, g5, g6, g7]

    # set up some 'reads'
    e1 = GenomicInterval("chr1", 2, 3)
    e2 = GenomicInterval("chr1", 4, 5)
    e3 = GenomicInterval("chr1", 4, 5)
    e4 = GenomicInterval("chr1", 7, 9)
    e5 = GenomicInterval("chr1", 7, 8)
    e6 = GenomicInterval("chr1", 10, 14)
    e7 = GenomicInterval("chr1", 12, 20)
    e8 = GenomicInterval("chr1", 20, 23)
    e9 = GenomicInterval("chr1", 24, 26)
    e10 = GenomicInterval("chr2", 1, 3)
    e11 = GenomicInterval("chr2", 2, 6)
    e12 = GenomicInterval("chr2", 4, 8)
    e13 = GenomicInterval("chr2", 10, 15)
    elements = [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13]

    # expected result
    expect = [(g1, set([e1, e2, e3])),
              (g2, set([e2, e3, e4, e5])),
              (g3, set([e4, e5, e6, e7])),
              (g4, set([e6, e7])),
              (g5, set([e7, e8])),
              (g6, set([e10, e11, e12, e13])),
              (g7, set([e11, e12]))]

    actual = [(x, set(l)) for x, l in bucketIterator(elements, buckets)]
    # compare against the interval tree approach too, just for good measure.
    trees = intervalTreesFromList(elements, openEnded=True)
    treeRes = [(g, set(trees[g.chrom].intersectingInterval(g.start, g.end)))
               for g in buckets]

    if debug:
      actualStr = "\n".join([str(x) + " == " + ",".join([str(y) for y in l])
                             for x, l in actual])
      treeReStr = "\n".join([str(x) + " == " + ",".join([str(y) for y in l])
                             for x, l in treeRes])
      expectStr = "\n".join([str(x) + " == " + ",".join([str(y) for y in l])
                             for x, l in expect])

      sys.stderr.write("actual result\n")
      sys.stderr.write(actualStr + "\n")
      sys.stderr.write("------------------------\n")
      sys.stderr.write("interval tree result\n")
      sys.stderr.write(treeReStr + "\n")
      sys.stderr.write("------------------------\n")
      sys.stderr.write("expected result\n")
      sys.stderr.write(expectStr + "\n")

    self.assertEqual(actual, expect)
    self.assertEqual(actual, treeRes)

  def testDistance(self):
    """
    Test getting absolute distance between intervals.

    :TODO: test needs to be implemented. Should test same and diff chroms.
    """
    pass

  def testSignedDistance(self):
    """
    Test getting signed distance between intervals.

    :TODO: test needs to be implemented. Should test same and diff chroms.
    """
    pass

  def testIntersects(self):
    """Test predicate for region overlap."""
    # region A contained entirely inside region B
    a = GenomicInterval("chr1", 10, 20)
    b = GenomicInterval("chr1", 15, 18)
    self.assertEqual(a.intersects(b), True)

    # region B contained entirely inside region A
    a = GenomicInterval("chr1", 15, 18)
    b = GenomicInterval("chr1", 10, 20)
    self.assertEqual(a.intersects(b), True)

    # region A overlaps start of region B
    a = GenomicInterval("chr1", 15, 20)
    b = GenomicInterval("chr1", 18, 25)
    self.assertEqual(a.intersects(b), True)

    # region A overlaps end of region B
    a = GenomicInterval("chr1", 18, 25)
    b = GenomicInterval("chr1", 10, 20)
    self.assertEqual(a.intersects(b), True)

    # region B overlaps one posiiton in A -- start
    a = GenomicInterval("chr1", 15, 18)
    b = GenomicInterval("chr1", 10, 16)
    self.assertEqual(a.intersects(b), True)

    # region B overlaps one position in A -- end
    a = GenomicInterval("chr1", 15, 18)
    b = GenomicInterval("chr1", 17, 26)
    self.assertEqual(a.intersects(b), True)

    # regions overlap on start/end, but different chroms
    a = GenomicInterval("chr1", 15, 18)
    b = GenomicInterval("chr2", 10, 20)
    self.assertEqual(a.intersects(b), False)

    # boundary case: region A ends immediately before B starts
    a = GenomicInterval("chr1", 15, 18)
    b = GenomicInterval("chr1", 18, 25)
    self.assertEqual(a.intersects(b), False)

    # boundary case: region A start immediately after B ends
    a = GenomicInterval("chr1", 16, 18)
    b = GenomicInterval("chr1", 10, 16)
    self.assertEqual(a.intersects(b), False)

    # boundary case: region B ends immediately before A starts
    a = GenomicInterval("chr1", 15, 18)
    b = GenomicInterval("chr1", 10, 15)
    self.assertEqual(a.intersects(b), False)

    # boundary case: region B start immediately after A ends
    a = GenomicInterval("chr1", 15, 18)
    b = GenomicInterval("chr1", 18, 25)
    self.assertEqual(a.intersects(b), False)

  def test_transform_center(self):
    """Test growing and shrinking regions around their center."""
    a = GenomicInterval("chr1", 15, 30, strand="+")
    b = GenomicInterval("chr1", 15, 30, strand="-")

    # symmetrical grow/shrink on positive and negative strand are equivalent
    # 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
    #                      --
    a.transform_center(11)
    b.transform_center(11)
    self.assertEqual(a.start, b.start)
    self.assertEqual(a.start, 17)
    self.assertEqual(a.end, b.end)
    self.assertEqual(a.end, 28)
    a.transform_center(15)
    b.transform_center(15)
    self.assertEqual(a.start, b.start)
    self.assertEqual(a.start, 15)
    self.assertEqual(a.end, b.end)
    self.assertEqual(a.end, 30)

    # non-symetric grow/shrink on positive and negative strand skewed, but
    # reverse gets you back to where you started.
    a.transform_center(10)
    b.transform_center(10)
    self.assertEqual((a.start, a.end), (17, 27))
    self.assertEqual((b.start, b.end), (18, 28))
    a.transform_center(15)
    b.transform_center(15)
    self.assertEqual((a.start, a.end), (15, 30))
    self.assertEqual((b.start, b.end), (15, 30))


###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
