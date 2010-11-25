#!/usr/bin/python

from Queue import Queue
import sys, unittest

class GenomicMapError(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)

class GenomicMap:
  """
    @summary: A genomic map translates chromosome, start, end tuples to
              a value. Entries without an end value are assume to end
              a single position after the start
  """
  def __init__(self):
    """
      @summary: GenomicMap constructor -- actual data is stored in vals, ages 
                allows us the keep track of the age of items (for flushing 
                later if they get too old)
    """
    self.vals = {}
    self.size = 0
    self.ages = Queue()
      
  def addValue(self, val, chrom, start, end = None):
    """
      @summary: add a value to the GenomicMap
      @param val: The value to add (can be anything)
      @param chrom: The chromosome that the value is associated with
      @param start: The start position that the value is associated with
      @param end: The end position that the value is associated with -- if
                  this is None then we assume an end = start + 1
    """
    if end == None : end = start + 1
    if end < start : 
      raise GenomicMapError("Trying to add element to GenomicMap which " +\
                            "begins (" + str(start) + ") after it ends " +\
                            "(" + str(end) + ")")
    if not chrom in self.vals : self.vals[chrom] = {}
    if (start,end) in self.vals[chrom] :
      raise GenomicMapError("GenomicMap already contains an entry for " +\
                            chrom + " from " + str(start) + " to " + str(end))
    self.vals[chrom][(start,end)] = val
    self.ages.put((chrom, start, end))
    self.size += 1
    
  def flush(self, amount = None):
    """
      @summary: remove the oldest <amount> items. If <amount> is None, or
                <amount> >= number of items in the map then remove all items
    """
    if amount == None or len(self) <= amount :
      self.vals = {}
      self.size = 0
      self.ages = Queue()
    else :
      removed = 0
      while removed < amount :
        if self.ages.empty() : break
        try :
          chrom, start, end = self.ages.get() 
          self.removeValue(chrom, start, end)
          removed += 1
        except GenomicMapError :
          pass
    
  def getValue(self, chrom, start, end = None):
    if end == None : end = start + 1
    if (not chrom in self.vals) or (not (start,end) in self.vals[chrom]) :
      raise GenomicMapError("no entry for " + chrom + " from " + str(start) +\
                            " to " + str(end))
    return self.vals[chrom][(start,end)]
  
  def removeValue(self, chrom, start, end = None):
    if end == None : end = start + 1
    if (not chrom in self.vals) or (not (start,end) in self.vals[chrom]) :
      raise GenomicMapError("no entry for " + chrom + " from " + str(start) +\
                            " to " + str(end))
    del self.vals[chrom][(start,end)]
    self.size -= 1
  
  def size(self, chrom = None):
    if chrom == None : return len(self)
    return len(self.vals[chrom])
  
  def __len__(self):
    return self.size
  
  def __str__(self):
    res = ""
    for chrom in self.vals :
      for start,end in self.vals[chrom] :
        v = self.vals[chrom][(start,end)]
        res += (chrom + ", " + str(start) + " -- " + str(end) + " --> " + str(v) + "\n")
    return res
  
class GetGenesTests(unittest.TestCase):
  """
    Unit tests for get genes 
  """
    
  def setUp(self):
    pass
  
  def testGenomicMap(self):
    debug = False
    
    m = GenomicMap()
    m.addValue("stuff", "chr1", 1, 5)
    m.addValue("stuff", "chr1", 7, 8)
    m.addValue("stuff", "chr1", 10, 15)
    m.addValue("stuff", "chr1", 16, 22)
    m.addValue("stuff", "chr2", 25, 28)
    
    self.assertTrue(len(m) == 5)
    self.assertRaises(GenomicMapError, m.addValue, "whatever", "chr1", 1, 5)
    self.assertRaises(GenomicMapError, m.removeValue, "chr1", 30, 40)
    m.removeValue("chr1", 10, 15)
    self.assertTrue(len(m) == 4)
    if debug :
      sys.stderr.write("contents before flush:\n")
      sys.stderr.write(str(m) +"\n")
    m.flush(3)
    if debug :
      sys.stderr.write("contents after flush:\n")
      sys.stderr.write(str(m) +"\n")
    self.assertTrue(len(m) == 1)
    self.assertTrue(str(m).strip() == "chr2, 25 -- 28 --> stuff")
  
if __name__ == "__main__":
    unittest.main(argv = [sys.argv[0]])