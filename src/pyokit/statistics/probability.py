#!/usr/bin/python

"""
  Date of Creation: 24th November 2010
  Description:      Classes and functions for generating and selecting
                    items based on probabilities

  Copyright (C) 2010-2014
  Philip J. Uren

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

import random, unittest, sys
from datastruct.intervalTree import IntervalTree

class WeightedRandom:
  """
    @summary Given a list of objects and a weight for each object, this class
             allows the random choice of an object (with replacement) such that
             the propbability of getting an object is proportional to it weight
             relative to the other objects.
  """

  class Interval :
    def __init__(self, start, end, obj):
      self.start = start
      self.end = end
      self.obj = obj
    def __str__(self):
      return str(self.start) + " to " + str(self.end) + " = " + str(self.obj)

  def __init__(self, candidates, weights):
    if len(candidates) != len(weights) :
      raise ProbabilityError("number of weights doesn't equal number of " +\
                             "objects")
    self._intervalTree = self._buildTree(weights, candidates)
    self._maxVal = sum(weights)

  def _buildTree(self, weights, candidates):
    """
      @summary: build interval tree from cumulative weights
    """
    intervals = []
    total = 0.0
    for i in range(0, len(weights)) :
      weight = weights[i]
      obj = candidates[i]
      start = total
      end = total + weight
      intervals.append(WeightedRandom.Interval(start, end, obj))
      total = total + weight
    return IntervalTree(intervals)

  def choose(self):
    rnd = random.uniform(0, self._maxVal)
    objs = self._intervalTree.intersectingPoint(rnd)
    if len(objs) != 1 :
      raise ProbabilityError("random choice failed, " + str(len(objs)) +\
                             " objects in same cumulative probability range!")
    return objs[0].obj


def generateProbabilities(num):
  """
    @summary: generate <num> probabilities such that for each probability x
              0 <= x <= 1 and sum(xs) = 1
    @param num: the number of probabilities to generate (int)
  """
  rnds = [random.random() for i in range(0,num)]
  rnds.sort()
  vals = []
  p = 0
  for n in rnds :
    vals.append(n - p)
    p = n
  vals.append(1 - p)
  return vals


class ProbabilityTests(unittest.TestCase):
  """
    @summary: Unit tests for probability
  """

  def setUp(self):
    pass

  def testWeightedChoice(self):
    numIter = 100000
    precision = 2

    items = ["a","b","c"]
    weights = [0.1, 0.3, 0.6]
    counts = {"a":0, "b":0, "c":0}
    r = WeightedRandom(items, weights)
    for i in range(numIter) :
      counts[r.choose()] += 1

    self.assertAlmostEqual(counts["a"] / float(numIter), 0.1, precision)
    self.assertAlmostEqual(counts["b"] / float(numIter), 0.3, precision)
    self.assertAlmostEqual(counts["c"] / float(numIter), 0.6, precision)

  def testGenerateProbabilities(self):
    numIter = 1000
    maxProbs = 100
    for i in range(0, numIter) :
      numProbs = int(random.random() * maxProbs)
      probs = generateProbabilities(numProbs)
      self.assertTrue(sum(probs) == 1)
      for prob in probs :
        self.assertTrue(prob >= 0)
        self.assertTrue(prob <= 1)


if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])
