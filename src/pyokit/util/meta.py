"""
  Date of Creation: 23rd April 2015

  Description:   General purpose utilities for meta programming

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
import inspect
import unittest


###############################################################################
#                             EXCEPTION CLASSES                               #
###############################################################################

class MetaError(Exception):
  """
  ...
  """
  def __init__(self, msg):
    self.value = msg

  def __str__(self):
    return repr(self.value)


###############################################################################
#                                  DECORATORS                                 #
###############################################################################

def decorate_all_methods(decorator):
  """
  This is a class decorator function that will apply the passed decorator to
  all of the methods in the decorated class, except the __init__ method. Its
  not specifically related to file indexing, but since that is the only place
  where I've used it and I have no better place to keep it at the moment...
  """
  def decorate_class(cls):
    for name, m in inspect.getmembers(cls, inspect.ismethod):
      if name != "__init__":
        setattr(cls, name, decorator(m))
    return cls
  return decorate_class


def decorate_all_methods_and_properties(method_decorator, property_decorator):
  """
  ...
  """
  def decorate_class(cls):
    for name, m in inspect.getmembers(cls, inspect.ismethod):
      if name != "__init__":
        setattr(cls, name, method_decorator(m))
    for name, p in inspect.getmembers(cls, lambda x: isinstance(x, property)):
      if name != "__init__":
        setattr(cls, name, property_decorator(name, p))
    return cls
  return decorate_class


def just_in_time_method(func):
  """
  This is a dcorator for methods. It redirect calls to the decorated method
  to the equivalent method in a class member called 'item'. 'item' is expected
  to be None when the class is instantiated. The point is to easily allow
  on-demand construction or loading of large, expensive objects just-in-time.

  To apply this decorator to a method in a class, the class must have the
  following instance variables:

  +--------------------------------------------------------------------------+
  | Name    |  Description                                                   |
  +=========+================================================================+
  | item    | The wrapped object. Calls to this method will be redirected    |
  |         | to the method with the same name in item. This should be set   |
  |         | to None when the object is created; it will be loaded          |
  |         | on-demand the first time this method is called.                |
  +---------+----------------------------------------------------------------+
  | factory | An object from which 'item' can be loaded. Can be a factory,   |
  |         | or similar, but must provide the subscript operator, as this   |
  |         | is used to pass a key that uniquely identifies 'item'          |
  +---------+----------------------------------------------------------------+
  | key     | any object that uniquely identifies 'item'; must be what is    |
  |         | expected by 'index' as argument for the subscript operator.    |
  +--------------------------------------------------------------------------+

  A common pattern is to create a new class B which inherits from A, implement
  the above requirements in B and then apply this decorator to all the methods
  inherited from A. If 'item' is an object of type A, then this pattern makes
  B behave exactly like A, but with just-in-time construction.
  """

  if not inspect.ismethod:
    raise MetaError("oops")

  def wrapper(self, *args, **kwargs):
    if self.item is None:
      self.item = self.factory[self.key]
    return getattr(self.item, func.__name__)(*args, **kwargs)
  return wrapper


def just_in_time_property(name, prop):
  """
  """
  if not isinstance(prop, property):
    raise MetaError("oops")

  def fget(self):
    if self.item is None:
      self.item = self.factory[self.key]
    return getattr(self.item, name)

  return property(fget)







###############################################################################
#                                 UNIT TESTS                                  #
###############################################################################
class TestMeta(unittest.TestCase):

  def test_just_in_time_construction(self):
    """
    Test the just-in-time decorator to wrap all methods in a class
    """
    class Record(object):
      def __init__(self, id, letters):
        self.id = id
        self.letters = letters

      def do_stuff(self, a):
        return a + len(self.letters)

      def __str__(self):
        return str(self.id) + " " + ",".join(self.letters)

    mapping = {1: Record(1, ["A", "B", "C", "D", "E"]),
               2: Record(2, ["F", "G", "H", "I", "J"]),
               3: Record(3, ["K", "L", "M", "N", "O"]),
               4: Record(4, ["P", "Q", "R", "S", "T"]),
               5: Record(5, ["U", "V", "W", "X", "Y"]),
               6: Record(6, ["Z", "A", "B", "C", "D"]),
               7: Record(7, ["E", "F", "G", "H", "I"]),
               8: Record(8, ["J", "K", "L", "M", "N"]),
               9: Record(9, ["O", "P", "Q", "R", "S"]),
               10: Record(10, ["T", "U", "V", "W", "X"])}

    class Factory(object):
      def __init__(self):
        pass

      def __getitem__(self, k):
        return mapping[k]

    def dummy_hash_2(item):
      return item.id

    @decorate_all_methods(just_in_time_method)
    class B(Record):
      def __init__(self, factory, key):
        self.factory = factory
        self.key = key
        self.item = None

    b1 = B(mapping, 5)
    b2 = B(mapping, 1)

    self.assertEqual(b1.do_stuff(12), 17)
    self.assertEqual(b2.do_stuff(12), 17)
    self.assertEqual(str(b1), "5 U,V,W,X,Y")  # U V W X Y
    self.assertEqual(str(b2), "1 A,B,C,D,E")  # A B C D E
