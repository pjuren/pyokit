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
import os
import sys

# pyokit imports
from pyokit.util.progressIndicator import ProgressIndicator


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class IndexError(Exception):

  """Exceptions that are raised when manipulating indexed files."""

  def __init__(self, msg):
    """:param msg: message for this exception."""
    self.value = msg

  def __str__(self):
    """:return: a string representation of this exception."""
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
    self._no_reindex = False
    try:
      # try treating this as a filename
      self.indexed_file = (indexed_file, None)
    except TypeError:
      try:
        # try treating this as a file handle
        self.indexed_file = (None, indexed_file)
      except TypeError:
        raise IndexError("failed to create index for " + str(indexed_file) +
                         "; reason: expected filename or stream-like "
                         "object, got " + str(type(indexed_file)))

  @property
  def indexed_file(self):
    """
    Getter for information on the file that this object indexes.

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
       (self.record_iterator is None or self.record_hash_function is None)):
      raise IndexError("Setting index file failed; reason: iterator "
                       "(self.record_iterator) or hash function "
                       "(self.record_hash_function) have to be set first")
    self._indexed_filename = filename
    self._indexed_file_handle = handle

  def __build_index(self, until=None, flush=False, verbose=False):
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

    if verbose:
      self._indexed_file_handle.seek(0, 2)  # seek to end
      total = self._indexed_file_handle.tell() - file_loc
      self._indexed_file_handle.seek(file_loc)  # back to where we were
      pind = ProgressIndicator(totalToDo=total,
                               messagePrefix="completed",
                               messageSuffix="of building out index")

    for item in self.record_iterator(self._indexed_file_handle):
      hash_val = self.record_hash_function(item)
      self._index[hash_val] = file_loc
      file_loc = self._indexed_file_handle.tell()
      if until is not None and hash_val == until:
        break
      if verbose:
        pind.done = file_loc
        pind.showProgress()

  def __len__(self):
    """:return: the current size (number of keys) in the index."""
    return len(self._index)

  def __getitem__(self, hash_value):
    """
    retrieve the record from the file which hashed to the given value.

    :return: the item in the indexed file that hashed to the given value.
    :raise IndexError:
    """
    # if we haven't seen this one yet, expand the index until we do..
    if hash_value not in self._index:
      if self._no_reindex:
        raise IndexError("No such hash value in index: " + str(hash_value))
      self.__build_index(until=hash_value)

    # if we still haven't seen it, fail; it cannot be retrieved (at least,
    # not at the moment)
    if hash_value not in self._index:
      fn = (" (" + str(self._indexed_filename) + ") "
            if self._indexed_filename is not None else "")
      raise IndexError("Unable to retrieve record (" + str(hash_value) +
                       ") from indexed file" + fn +
                       "; reason: no such record found")

    # we got a match to the hash_value, attempt to seek to this location and
    # return the next element
    orig_pos = self._indexed_file_handle.tell()
    self._indexed_file_handle.seek(self._index[hash_value])
    r_iter = self.record_iterator(self._indexed_file_handle)
    try:
      return r_iter.next()
    except StopIteration:
      fn = (" (for " + str(self._indexed_filename) + ")"
            if self._indexed_filename is None else "")
      raise IndexError("Fatal error, index" + fn +
                       " specifies location beyond end of indexed file")
    finally:
      self._indexed_file_handle.seek(orig_pos)

  def __eq__(self, other):
    """
    determine whether two IndexedFile objects are equal. To be so, the index
    dictionaries themselves must be equal. We will not compare filename;
    basically, we don't care where the index information came from or whether
    it might cease to be equal later, as long as the actual index itself is
    the same.

    :param other: the other IndexedFile to compare against.
    """
    return self._index == other._index

  def __str__(self):
    """
    Produce a string representation of this index.

    Use this with caution, as the full index is converted to string, which
    might be quite large.
    """
    res = ""
    for key in self._index:
      res += (str(key) + "\t" + str(self._index[key]) + "\n")
    return res

  def write_index(self, fh, to_str_func=str, generate=True, verbose=False):
    """
    Write this index to a file.

    Only the index dictionary itself is stored, no informatiom about the
    indexed file, or the open filehandle is retained. The Output format is
    just a tab-separated file, one record per line. The last column is the
    file location for the record and all columns before that are collectively
    considered to be the hash key for that record (which is probably only 1
    column, but this allows us to permit tabs in hash keys).

    :param fh:           either a string filename or a stream-like object to
                         write to.
    :param to_str_func:  a function to convert hash values to strings. We'll
                         just use str() if this isn't provided.
    :param generate:     build the full index from the indexed file if it
                         hasn't already been built. This is the default, and
                         almost certainly what you want, otherwise just the
                         part of the index already constructed is written
                         (which might be nothing...)
    :param verbose:      if True, output progress messages to stderr.
    """
    try:
      handle = open(fh, "w")
    except TypeError:
      # okay, not a filename, try to treat it as a stream to write to.
      handle = fh
    if generate:
      self.__build_index(verbose=verbose)
    for key in self._index:
      handle.write(to_str_func(key) + "\t" + str(self._index[key]) + "\n")

  def __iter__(self):
    """:return: an iterator over the index keys."""
    return self._index.__iter__()

  def read_index(self, fh, indexed_fh, rec_iterator=None,
                 rec_hash_func=None, parse_hash=str, flush=True,
                 no_reindex=True, verbose=False):
    """
    Populate this index from a file. Input format is just a tab-separated file,
    one record per line. The last column is the file location for the record
    and all columns before that are collectively considered to be the hash key
    for that record (which is probably only 1 column, but this allows us to
    permit tabs in hash keys). Lines consisting only of whitespace are skipped.

    :param fh:            filename or stream-like object to read from.
    :param indexed_fh:    either the filename of the indexed file or handle to
                          it.
    :param rec_iterator:  a function that will return an interator for the
                          indexed file type (not the iterator for the file
                          itself). This function must take a single argument
                          which is the name the file to iterate over, or a
                          stream like object similar to a filestream.
    :param rec_hash_func: a function that accepts the record type produced by
                          the iterator and produces a unique hash for each
                          record.
    :param parse_hash:    a function to convert the string representation of
                          the hash into whatever type is needed. By default,
                          we just leave these as strings.
    :param flush:         remove everything currently in the index and discard
                          any details about a file that is already
                          fully/partially indexed by this object. This is the
                          default behavior. If False, then data from <fh> is
                          just added to the existing index data (potentially
                          overwriting some of it) and the existing index can
                          continue to be used as before.
    :param no_reindex:    if True, after loading the index, a missing key will
                          cause an exception, rather than trigger re-scanning
                          the indexed file for the associated record. The only
                          reason to set this to False would be if your index
                          was incomplete.
    :param verbose:       output status message to STDERR about progress
                          reading the index (if possible).

    :raise IndexError: on malformed line in input file/stream
    """
    # set the record iterator and hash functions, if they were given
    if rec_iterator is not None:
      self.record_iterator = rec_iterator
    if rec_hash_func is not None:
      self.record_hash_function = rec_hash_func

    # disable re-indexing?
    self._no_reindex = no_reindex

    # figure out what kind of index identifier we got: handle or filename?
    handle = fh
    try:
      handle = open(fh)
    except TypeError:
      # okay, not a filename, we'll try treating it as a stream to read from.
      pass

    # clear this index?
    if flush:
      self._index = {}
      self._indexed_file_handle = None
      self._indexed_file_name = None

    # replace the name/handle for the indexed file
    indexed_fn = None
    try:
      # try treating this as a filename
      self.indexed_file = (indexed_fh, None)
      indexed_fn = indexed_fh
    except TypeError:
      try:
        # try treating this as a file handle
        self.indexed_file = (None, indexed_fh)
      except TypeError:
        fn = " from " + str(fh) if indexed_fn is not None else ""
        raise IndexError("failed to read index" + fn + "; "
                         "reason: expected indexed filename or stream-like "
                         "object, got " + str(type(indexed_fh)))

    # try to get an idea of how much data we have...
    if verbose:
      try:
        total = os.path.getsize(handle.name)
        pind = ProgressIndicator(totalToDo=total, messagePrefix="completed",
                                 messageSuffix="of loading " + handle.name)
      except AttributeError as e:
        sys.stderr.write(str(e))
        sys.stderr.write("completed [unknown] of loading index")
        verbose = False

    # read the index file and populate this object
    for line in handle:
      line = line.rstrip()

      if verbose:
        pind.done = handle.tell()
        pind.showProgress()

      if line.isspace():
        continue
      parts = line.split("\t")
      if len(parts) < 2:
        raise IndexError("failed to parse line: '" + line + "'")
      key = parse_hash("\t".join(parts[:-1]))
      value = parts[-1]
      self._index[key] = int(value)


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

class TestIndexedFile(unittest.TestCase):

  def setUp(self):
    """
    Test indexing a simple file format.

    Records in the indexed files have the following structure: 3 A B C D E
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
    """Test grabbing an item when nothing has been indexed yet."""
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
    """Test grabbing an item that has already been indexed."""
    index = IndexedFile(StringIO.StringIO(self.test_case_0), self.r_iter,
                        self.r_hash)
    self.assertEqual(index[5], (5, ["U", "V", "W", "X", "Y"]))
    self.assertEqual(index[5], (5, ["U", "V", "W", "X", "Y"]))

  def test_indexedFile_not_present(self):
    """Test asking for an item that does not exist in the file."""
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
    """Test writing an index to file, reading it back and checking equality."""
    # create dummy files
    indexed_fh = StringIO.StringIO(self.test_case_0)
    index_fh = StringIO.StringIO()

    # index the dummy 'indexed' file to create an in-memory index
    index = IndexedFile(indexed_fh, self.r_iter, self.r_hash)

    # write out in-memory index to 'disk'
    index.write_index(index_fh)

    # reset both indexed and index files back to start, then read index
    # from 'disk' into new index object
    index_fh.seek(0)
    index_2 = IndexedFile(record_iterator=self.r_iter,
                          record_hash_function=self.r_hash)
    index_2.read_index(index_fh, indexed_fh, parse_hash=int)

    # compare...
    self.assertEqual(index[5], index_2[5])
    self.assertEqual(index, index_2)


###############################################################################
#                 ENTRY POINT IF RUN AS STAND-ALONE MODULE                    #
###############################################################################

if __name__ == '__main__':
    unittest.main()
