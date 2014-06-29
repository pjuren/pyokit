#!/usr/bin/python

"""
Date of Creation: 3rd June 2010
An Interval Tree is a data structure for quickly determining the set of
intervals that intersect a given point or interval.

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

import copy, random, unittest

class IntervalTreeError(Exception):
  """
    Exception class for errors that occur when building or using interval trees
  """
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)


class IntervalTreeNode :
  def __init__(self, intervals, mid):
    self.starts = copy.copy(intervals)  # shallow copies
    self.ends = copy.copy(intervals)
    self.starts = sorted(self.starts, key=lambda x : x.start)
    self.ends = sorted(self.ends, key=lambda x : x.end)
    self.mid = mid

  def __str__(self):
    res = ""
    for item in self.starts :
      res += "(" + str(item.start) + "--" + str(item.end) + ")"
    return res


class IntervalTree :
  """
    An interval tree is a binary tree that allows fast O(log(n)) lookup of
    intervals that intersect a given point or interval.

    :param intervals: list of intervals, doesn't need to be sorted in any way.
                      Can be any object, as long as they have 'start' and 'end'
                      attributes.
  """
  def __init__(self, intervals, openEnded=False):
    self.openEnded = openEnded
    self.left = None
    self.right = None

    ## can't build a tree with no intervals...
    if intervals == None or len(intervals) <= 0:
      raise IntervalTreeError("Interval tree constructor got empty " +\
                              "or null set of intervals")

    # sort the list by start index (this is arbitrary and just lets
    # us try to split it evenly)
    intervals = sorted(intervals, key=lambda x : x.start)

    # pick a mid-point and split the list
    mid = intervals[len(intervals) / 2]
    mid = int((mid.end - mid.start) / 2) + mid.start

    here = []
    lt = []
    rt = []
    for i in intervals :
      # place all intervals that end before <mid> into the left subtree
      if (not self.openEnded and i.end < mid) or \
         (self.openEnded and i.end <= mid) :
        lt.append(i)
      # place all intervals that begin after <mid> into the right subtree
      elif i.start > mid :
        rt.append(i)
      # everything else must overlap mid, so we keep it here
      else :
        here.append(i)

    if len(here) <= 0 :
      intStrs = ",".join(str(x.start) + " -- " + str(x.end) for x in intervals)
      raise IntervalTreeError("picked mid point at " + str(mid) +\
                              " but failed to intersect any of " + intStrs)

    if len(lt) > 0 : self.left = IntervalTree(lt, self.openEnded)
    if len(rt) > 0 : self.right = IntervalTree(rt, self.openEnded)
    self.data = IntervalTreeNode(here, mid)


  def intersectingPoint(self, p):
    """
    given a point, determine which set of intervals in the tree are intersected.

    :param p: intersection point
    :return: the list of intersected intervals
    """

    # perfect match
    if p == self.data.mid :
      return self.data.ends

    if p > self.data.mid :
      # we know all intervals in self.data begin before p (if they began after
      # p, they would have not included mid) we just need to find those that
      # end after p
      endAfterP = [r for r in self.data.ends \
                   if (r.end >= p and not self.openEnded) or \
                   (r.end > p and self.openEnded)]
      if self.right != None : endAfterP.extend(self.right.intersectingPoint(p))
      return endAfterP

    if p < self.data.mid :
      # we know all intervals in self.data end after p (if they ended before p,
      # they would have not included mid) we just need to find those that start
      # before p
      startBeforeP = [r for r in self.data.starts if r.start <= p]
      if self.left != None : startBeforeP.extend(self.left.intersectingPoint(p))
      return startBeforeP


  def intersectingInterval(self, start, end):
    """
      given an interval, determine which set of intervals in the tree are
      intersected.

      :param start: start of the intersecting interval
      :param end:   end of the intersecting interval
      :return:      the list of intersected intervals
    """

    # find all intervals in this node that intersect start and end
    l = []
    for x in self.data.starts :
      xStartsAfterInterval = (x.start > end and not self.openEnded) or \
                             (x.start >= end and self.openEnded)
      xEndsBeforeInterval = (x.end < start and not self.openEnded) or \
                            (x.end <= start and self.openEnded)
      if ((not xStartsAfterInterval) and (not xEndsBeforeInterval)) :
        l.append(x)

    # process left subtree (if we have one) if the requested interval begins
    # before mid
    if self.left != None and start <= self.data.mid:
      l += self.left.intersectingInterval(start, end)

    # process right subtree (if we have one) if the requested interval ends
    # after mid
    if self.right != None and end >= self.data.mid:
      l += self.right.intersectingInterval(start, end)

    return l

  def intersectingIntervalIterator(self, start, end):
    """
      Get an iterator which will iterate over those objects in the tree which
      intersect the given interval - sorted in order of start index

      :param start: find intervals in the tree that intersect an interval with
                    with this start index (inclusive)
      :param end:   find intervals in the tree that intersect an interval with
                    with this end index (exclusive)
      :return: an iterator that will yield intersected intervals
    """
    items = self.intersectingInterval(start, end)
    items.sort(key = lambda x: x.start)
    for item in items :
      yield item

  def __str__(self):
    return str(self.left) + "," + str(self.data) + "," + str(self.right)


class TestIntervalTree(unittest.TestCase):
  class TestInterval :
    """ a small internal class used only for testing """
    def __init__(self, s, e):
      self.start = s
      self.end = e
    def str(self):
      return str(self.start) + " - " + str(self.end)
    def __repr__(self):
      return str(self.start) + " - " + str(self.end)

  def setUp(self):
    self.testIntervals = []
    self.NUM_TEST_INTERVALS = 1000
    self.MAX_INTERVAL = 10000
    self.MIN_INTERVAL = 0
    self.NUM_TESTS = 100

    ## create some test intervals
    for i in range (0, self.NUM_TEST_INTERVALS) :
      s = (random.random() * (self.MAX_INTERVAL - self.MIN_INTERVAL)) +\
          self.MIN_INTERVAL
      e = (random.random() * (self.MAX_INTERVAL - s)) + s
      self.testIntervals.append(TestIntervalTree.TestInterval(s,e))
    self.tree = IntervalTree(self.testIntervals)

  def testEndPoints(self):
    interval = TestIntervalTree.TestInterval(5,10)
    tree = IntervalTree([interval])

    self.assertTrue(len(tree.intersectingPoint(4)) == 0)
    self.assertTrue(len(tree.intersectingPoint(5)) == 1)
    self.assertTrue(len(tree.intersectingPoint(6)) == 1)
    self.assertTrue(len(tree.intersectingPoint(9)) == 1)
    self.assertTrue(len(tree.intersectingPoint(10)) == 1)
    self.assertTrue(len(tree.intersectingPoint(11)) == 0)

  def testEmpty(self):
    """ Test that the interval tree raises an exception if provided with an
    emtpy list
    """
    self.assertRaises(IntervalTreeError, IntervalTree, [])
    self.assertRaises(IntervalTreeError, IntervalTree, None)


  def testIntersectingPoint(self):
    for i in range(0,self.NUM_TESTS) :
      point = random.random() * (self.MAX_INTERVAL - self.MIN_INTERVAL) +\
              self.MIN_INTERVAL

      ## do it the slow way...
      correctAnswer = []
      for interval in self.testIntervals :
        if point > interval.start and point < interval.end :
          correctAnswer.append(interval)

      ## now use the interval tree...
      actualAnswer = self.tree.intersectingPoint(point)

      ## compare...
      correctAnswer.sort(key = lambda x: x.start)
      actualAnswer.sort(key = lambda x: x.start)
      self.assertEquals(correctAnswer, actualAnswer)


  def testIntersectingInterval(self):
    for i in range(0,self.NUM_TESTS) :
      # create a random interval..
      s = (random.random() * (self.MAX_INTERVAL - self.MIN_INTERVAL)) +\
          self.MIN_INTERVAL
      e = (random.random() * (self.MAX_INTERVAL - s)) + s

      ## work out the answer the slow way - look at every interval
      correctAnswer = []
      for interval in self.testIntervals :
        if (s > interval.start and s < interval.end) \
        or (e > interval.start and e < interval.end) \
        or (interval.start > s and interval.start < e) \
        or (interval.end > s and interval.end < e) :
          correctAnswer.append(interval)

      ## now use the tree..
      actualAnswer = self.tree.intersectingInterval(s, e)

      ## compare...
      correctAnswer.sort(key = lambda x: x.start)
      actualAnswer.sort(key = lambda x: x.start)
      self.assertEquals(correctAnswer, actualAnswer)

  def testThree(self):
    class testInterval:
      def __init__(self,s,e):
        self.start = s
        self.end = e
      def __eq__(self, other):
        return self.end == other.end and self.start == other.start
    dead1 = testInterval(215844, 545297)
    tree = IntervalTree([dead1])

    ans = tree.intersectingInterval(296532, 592921)
    self.assertEqual(ans, [dead1])

  def testOpenVsClosedInterval(self):
    """ test that passing the open interval switch results in intersections
    with the end of an interval being dropped
    """
    one = TestIntervalTree.TestInterval(10,15)
    two = TestIntervalTree.TestInterval(15,20)
    t = IntervalTree([one,two])
    t2 = IntervalTree([one,two], openEnded=True)

    self.assertTrue(len(t.intersectingPoint(15)) == 2)
    self.assertTrue(len(t2.intersectingPoint(15)) == 1)

if __name__ == '__main__':
    unittest.main()
