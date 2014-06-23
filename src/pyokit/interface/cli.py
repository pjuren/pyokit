#!/usr/bin/python

"""
  Date of Creation: 20th April 2010
  Description:      Command line parsing module, extends getopt functionality

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

import sys, copy, getopt

DEFAULT_TYPE = str

class InterfaceException(Exception):
  def __init__(self, msg):
    self.value = msg
  def __str__(self):
    return repr(self.value)


class Option :
  def __init__(self, short, long, description, argName = None,
               default = None, required = False, type = DEFAULT_TYPE,
               special = False):
    """
      @summary: Constructor for Option class
      @param short: the short name for the option (should be 1 char long)
      @param long: the long name for the option (arbitrary length)
      @param description: text description of what the option does
      @param argName: name of the argument, if the option requires one
      @param default: if the option is not given, what default value is used?
      @param required: will passing the cmd line fail if this option is missing?
      @param type: data type for this option (e.g. int, str, float, etc.)
      @param special: if True, disregard requirements for other options when
                      this option is set
    """
    self.short = short
    self.long = long
    self.description = description

    self.required = required
    self.argName = argName
    self.type = type
    self.special = special

    self.default = default
    self.value = None
    self.set = False

    # if a default value is given, we set it to that on creation
    # (note that we don't count this as 'setting' the option)
    if self.default != None :
      self.value = self.type(default)

    # short option names should only be 1 char
    if len(self.short) != 1 :
      raise InterfaceException("short option names must be 1 character long")

    # doens't make sense to tell us the type but not the name of the arg
    if self.type != DEFAULT_TYPE and argName == None :
      raise InterfaceException("Got argument type but no name")

  def isSet(self):
    return self.set

  def isRequired(self):
    return self.required

  def __str__(self):
    return  "-" + self.short + ", --" + self.long + " : " +\
            self.description.replace("\t","\\t") + " (required? " +\
            str(self.required) + ")"

class CLI :
  def __init__(self, progName, shortDesc, longDesc):
    self.programName = progName
    self.longDescription = longDesc
    self.shortDescription = shortDesc

    self.options = []

    self.minArgs = None
    self.maxArgs = None


  def usage(self):
    print self.programName + " [options]",
    if self.maxArgs == -1 or self.maxArgs - self.minArgs > 1 : print "[files..]"
    elif self.maxArgs - self.minArgs == 1: print "[file]"
    print

    print self.shortDescription

    print "Options:"
    for o in self.options :
      print "\t", ("-" + o.short + ", --" + o.long).ljust(15), ":",
      print o.description.replace("\t","\\t").ljust(80)

  def parseCommandLine(self, line):
    """
      @summary: Parse the given command line.
      @param line: list of tokens from the command line. Should have had
                   program name removed - generally one would do this:
                   sys.argv[1:]
    """

    # separate arguments and options
    try:
      optlist, args = getopt.getopt(line, self._optlist(), self._longoptl())
    except getopt.GetoptError, err:
      print str(err)
      self.usage()
      sys.exit(2)

    # parse options
    for opt, arg in optlist :
      opt = opt.lstrip("-")
      found = False
      for option in self.options :
        if option.short == opt or option.long == opt :
          option.value = option.type(arg)
          option.set = True
          found = True
      if not found :
        sys.stderr.write("unknown option: " + opt)
        self.usage()
        sys.exit(2)

    # check for special options
    specialOptions = False
    for opt in self.options :
      if opt.special and opt.isSet() : specialOptions = True

    # parse arguments
    if len(args) < self.minArgs and not specialOptions :
      err = str(len(args)) + " input arguments found, but required at least " +\
            str(self.minArgs) + " arguments"
      sys.stderr.write(err)
      sys.stderr.write("\n\n")
      self.usage()
      sys.exit(2)
    else :
      self.args = args

    # check for required options
    for option in self.options :
      if option.isRequired() and not option.isSet() and not specialOptions:
        err = "required option \'" + str(option) + "\' is missing"
        sys.stderr.write(err)
        print "\n"
        self.usage()
        sys.exit(2)

  def getOption(self, name):
    """
      @summary returns the Option associated with the given name - name can
               be short or long name.
      @raise InterfaceException: if the named option doesn't exist.
    """
    name = name.strip()
    for o in self.options :
      if o.short == name or o.long == name :
        return o
    raise InterfaceException("No such option: " + name)

  def getValue(self, name):
    """
      @summary: returns the value of the option matching the name given -
                can be long or short name
      @raise InterfaceException: if the named option doesn't exist.
    """
    return self.getOption(name).value

  def hasOption(self, name):
    """
      @summary: returns true if an option exists in the UI for the given
                name (short or long name ok)
    """
    name = name.strip()
    for o in self.options :
      if o.short == name or o.long == name :
        return True
    return False

  def addOption(self, o):
    """
      @summary: add a new Option to the UI
      @param o: the option to add
      @raise InterfaceException: if an option with that name already exists.
    """
    if self.hasOption(o.short) :
      raise InterfaceException("Failed adding option - already have " + str(o))
    self.options.append(o)

  def optionIsSet(self, name):
    """
      @summary: returns true if an option matching the given name exists and
                it has had it's value set by the user
    """
    name = name.strip()
    if not self.hasOption(name) : return False
    return self.getOption(name).isSet()

  def hasArgument(self, num):
    """
      @summary: Returns true if the user has supplied at least <num> + 1
                arguments. i.e. arguments are indexed from 0, so asking for
                the 0th arg will check to see if the user supplied 1 or
                more args
    """
    return len(self.args) >= num + 1

  def getArgument(self, num):
    """
      @summary: return the num^th argument. Arguments are always parsed as
                strings, so this will be a string. Indexing here is from 0.
      @param num: the argument number to fetch.
      @raise InterfaceException: if there are fewer than num + 1 arguments.
    """
    if not self.hasArgument(num) :
      raise InterfaceException("Failed to retrieve argument " + str(num) +\
                               " -- not enough arguments provided")
    return self.args[num]

  def getAllArguments(self):
    """
      @summary: return a list of all the arguments.
    """
    # don't leak a reference to the internal variable args
    return copy.copy(self.args)

  def _optlist(self):
    """
      @summary: TODO
    """
    res = ""
    for o in self.options :
        res += o.short
        if o.argName != None :
          res += ":"
    return res

  def _longoptl(self):
    """
      @summary: TODO
    """
    res = []
    for o in self.options :
      nm = o.long
      if o.argName is not None :
        nm += "="
      res.append(nm)
    return res
