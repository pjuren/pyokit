#!/usr/bin/python

""" 
  Date of Creation: 24th November 2010     
                      
  Description:        Classes and functions for generating and selecting
                      items based on probabilities

  Copyright (C) 2010  
  University of Southern California,
  Philip J. Uren,
  Andrew D. Smith
  
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
  
  --------------------
  
  Known Bugs:    None
  
  Revision 
  History:       None 
  
  
  TODO:          None 
"""

import random, unittest, sys

def weightedChoice(candidates, weights):
  """
    @summary: given a list of objects and a weight for each one, choose an item
              randomly with probability based on the weights
    @param candidates: list of objects to choose from
    @param weights: list of weights (float or int) for each item in candidates
                    weights need not be normalised 
  """
  cumWeights = []
  total = 0
  for weight in weights:
    total = total + weight
    cumWeights.append(total)
  rnd = random.uniform(0, total)
  for candidate, cumWeight in zip(candidates, cumWeights):
    if rnd <= cumWeight:
      return candidate

  
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
    for i in range(numIter) :
      counts[weightedChoice(items, weights)] += 1
     
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