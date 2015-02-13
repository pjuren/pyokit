"""
  Date of Creation: 12th Feb 2015

  Description:   Classes and related functions for indexing files that are made
                 up of many records, each of which can be treated independently
                 and can be represented by a unique identifier. Speeds up
                 random access to such files when storing the whole thing in
                 memory is not possible.

  Copyright (C) 2010-2015
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

# standard python imports
import StringIO
import unittest


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class IndexError(Exception):
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


###############################################################################
#                              INDEXING CLASSES                               #
###############################################################################

class IndexedFile(object):
  """
  An IndexedFile is a data structure that indexes a large file that is made up
  of many individual records. Any file for which an iterator exists to iterate
  over these records, and a hash function exists to extract a unique key from
  each record can be indexed.

  The reduction in size afforded by the index, as opposed to storing the full
  file in memory will be greatest when the records themselves are large and the
  hash values are small. If your application needs to access each record only
  once, then speed (after indexing) should be identical.

  The file is indexed on-demand; we never never index more than is needed to
  get the desired record. Of course, if the desireded record is not present,
  we will need to index the whole file to determine this.

  :param indexed_filename:     file that is being indexed; can be 'None', in
                               which case an 'empty' index is created. Can be
                               either a filename or a stream-like object that
                               behaves like a file handle.
  :param record_iterator:      a function that will return an interator for the
                               indexed file type (not the iterator for the
                               file itself). This function must take a single
                               argument which is the name the file to iterate
                               over, or a stream like object similar to a
                               filestream.
  :param record_hash_function: a function that accepts the record type produced
                               by the iterator and produces a unique hash for
                               each record.
  """

  def __init__(self, indexed_file=None, record_iterator=None,
               record_hash_function=None):
    self.record_iterator = record_iterator
    self.record_hash_function = record_hash_function
    self._index = {}
    self._indexed_filename = None
    self._indexed_file_handle = None
    try:
      # try treating this as a filename
      self.indexed_file = (indexed_file, None)
    except TypeError:
      try:
        # try treating this as a file handle
        self.indexed_file = (None, indexed_file)
      except TypeError:
        raise IndexError("failed to create index for " + str(indexed_file)
                         + "; reason: expected filename or stream-like "
                         "object, got " + str(type(indexed_file)))

  @property
  def indexed_file(self):
    """
    Getter for information on the file that this object indexes
    :return: the tuple (filename, handle) -- either of which might be None.
    """
    return (self._indexed_filename, self._indexed_file_handle)

  @indexed_file.setter
  def indexed_file(self, f):
    """
    Setter for information about the file this object indexes.
    :param f: a tuple of (filename, handle), either (or both) of which can be
              None. If the handle is None, but filename is provided, then
              handle is created from the filename. If both handle and filename
              are None, or they don't match the previous values indexed by this
              object, any current data in this index is cleared. If either are
              not None, we require the iterator and the hash function for this
              object to already be set.
    """
    filename, handle = f
    if handle is None and filename is not None:
      handle = open(filename)
    if (handle is None and filename is None) or \
       (filename != self._indexed_filename) or \
       (handle != self._indexed_file_handle):
      self.index = {}
    if ((handle is not None or filename is not None) and
       (self.record_iterator == None or self.record_hash_function == None)):
      raise IndexError("Setting index file failed; reason: iterator "
                       "(self.record_iterator) or hash function "
                       "(self.record_hash_function) have to be set first")
    self._indexed_filename = filename
    self._indexed_file_handle = handle

  def __build_index(self, until=None, flush=False):
    """
    build/expand the index for this file.
    :param until: expand the index until the record with this hash has been
                  incorporated and then stop. If None, go until the iterator
                  is exhausted. Note that if this hash is already in the index,
                  no new items will be
    :param flush: if True, anything already in the index is discarded.
    """
    assert(self._indexed_file_handle is not None)
    if flush:
      self._index = {}
    file_loc = self._indexed_file_handle.tell()
    for item in self.record_iterator(self._indexed_file_handle):
      hash_val = self.record_hash_function(item)
      self._index[hash_val] = file_loc
      file_loc = self._indexed_file_handle.tell()
      if until is not None and hash_val == until:
        break

  def __getitem__(self, hash_value):
    """
    retrieve the record from the file which hashed to the given value.
    :return: the item in the indexed file that hashed to the given value.
    :raise IndexError:
    """
    # if we haven't seen this one yet, expand the index until we do..
    if hash_value not in self._index:
      self.__build_index(until=hash_value)

    # if we still haven't seen it, fail; it cannot be retrieved (at least,
    # not at the moment)
    if hash_value not in self._index:
      fn = (" (" + str(self._indexed_filename) + ") "
            if self._indexed_filename is None else "")
      raise IndexError("Unable to retrieve record from indexed file" + fn
                       + "; reason: no such record found")

    # we got a match to the hash_value, attempt to seek to this location and
    # return the next element
    self._indexed_file_handle.seek(self._index[hash_value])
    r_iter = self.record_iterator(self._indexed_file_handle)
    try:
      return r_iter.next()
    except StopIteration:
      fn = (" (for " + str(self._indexed_filename) + ")"
            if self._indexed_filename is None else "")
      raise IndexError("Fatal error, index" + fn
                       + " specifies location beyond end of indexed file")

  def write_index(self, fn):
    """
    write this index to a file
    """
    raise IndexError("writing index to file not yet implemented")

  def read_index(self, fn):
    """
    populate this index from a file
    :raise IndexError: if the record iterator or record hash function is not
                       set
    """
    raise IndexError("reading index from file not yet implemented")


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestAlignmentIterators(unittest.TestCase):

  def setUp(self):
    """
    We will test indexing a simple file where records have the following
    structure:

      3 A B C D E

    Here, all data about each record is on a single line -- white-space
    separated. The first item on the data line is a unique ID for the element.
    """

    def dummy_iterator(strm):
      for line in strm:
        line = line.strip()
        if line == "":
          continue
        parts = line.split()
        yield (int(parts[0]), parts[1:])
    self.r_iter = dummy_iterator

    def dummy_hash(item_tuple):
      return int(item_tuple[0])
    self.r_hash = dummy_hash

    self.test_case_0 = "1  A B C D E\n" +\
                       "2  F G H I J\n" +\
                       "3  K L M N O\n" +\
                       "4  P Q R S T\n" +\
                       "5  U V W X Y\n" +\
                       "6  Z A B C D\n" +\
                       "7  E F G H I\n" +\
                       "8  J K L M N\n" +\
                       "9  O P Q R S\n" +\
                       "10 T U V W X"

  def test_indexedFile_raw(self):
    """
      test grabbing an item when nothing has been indexed yet.
    """
    index = IndexedFile(StringIO.StringIO(self.test_case_0), self.r_iter,
                        self.r_hash)
    self.assertEqual(index[5], (5, ["U", "V", "W", "X", "Y"]))
    index = IndexedFile(StringIO.StringIO(self.test_case_0), self.r_iter,
                        self.r_hash)
    self.assertEqual(index[1], (1, ["A", "B", "C", "D", "E"]))
    index = IndexedFile(StringIO.StringIO(self.test_case_0), self.r_iter,
                        self.r_hash)
    self.assertEqual(index[10], (10, ["T", "U", "V", "W", "X"]))

  def test_indexedFile_already_indexed(self):
    """
      test grabbing an item that has already been indexed.
    """
    index = IndexedFile(StringIO.StringIO(self.test_case_0), self.r_iter,
                        self.r_hash)
    self.assertEqual(index[5], (5, ["U", "V", "W", "X", "Y"]))
    self.assertEqual(index[5], (5, ["U", "V", "W", "X", "Y"]))

  def test_indexedFile_not_present(self):
    """
      test asking for an item that does not exist in the file
    """
    index = IndexedFile(StringIO.StringIO(self.test_case_0), self.r_iter,
                        self.r_hash)
    self.assertRaises(IndexError, index.__getitem__, 11)
    self.assertRaises(IndexError, index.__getitem__, 0)

  def test_indexedFile_allFullIndex(self):
    """
      test asking for an item after forcing a full index by first requesting
      an item that is not present.
    """
    index = IndexedFile(StringIO.StringIO(self.test_case_0), self.r_iter,
                        self.r_hash)
    self.assertRaises(IndexError, index.__getitem__, 11)
    self.assertEqual(index[6], (6, ["Z", "A", "B", "C", "D"]))

  def test_indexedFile_write_read_equality(self):
    """
      test writing an index to file, reading it back and checking equality.
      :note: these aren't implemented yet, when they are this test should be
             updated.
    """
    index = IndexedFile(StringIO.StringIO(self.test_case_0), self.r_iter,
                        self.r_hash)
    self.assertRaises(IndexError, index.write_index, StringIO.StringIO())
    self.assertRaises(IndexError, index.read_index, StringIO.StringIO())


###############################################################################
#                 ENTRY POINT IF RUN AS STAND-ALONE MODULE                    #
###############################################################################

if __name__ == '__main__':
    unittest.main()
