Pyokit: Python Bio Toolkit
==========================

Pyokit is a python library with a collection of tools to make the processing
and analysis of high-throughput biological datasets in python easier.

Dependencies
------------
Pyokit depends on rpy2 for some of its functionality. You don't strictly
need it, Pyokit can be installed without it, but if you don't have it some
functionality will not be supported.

Installing
----------
**From PyPI.** Pyokit is available from the Python package index. The easiest
way to install it is to use pip. Note that You may need to prefix certain
commands below with ```sudo``` if you're installing to the global python
installation.

1. Get pip (if you don't already have it):
```bash
  wget https://bootstrap.pypa.io/get-pip.py
  python get-pip.py
```
This will install pip for you. if ```python``` doesn't point to the version of
python you want to use, then replace it with the full path to the one you prefer
(you'll probably have to do that if you don't have admin. access on the machine
you're using).  
2. Now you can install pyokit:
```bash
pip install pyokit
```
... and you're done.

**From source.** The bleeding-edge version can always be cloned from the GitHub
repository (www.github.com/pjuren/pyokit), where tagged releases can also be
downloaded as .tar.gz. If you take this route, you will still need pip (see
above if you don't have it). After unpacking the distribution (or cloning the
repository), ```cd``` into the newly created directory and then
1. Build the distribution package:
```bash
python setup.py sdist
```
2. Install the distribution:
```bash
pip install --no-index dist/pyokit-a.b.c.tar.gz
```
where a.b.c is the version you downloaded.

Using Pyokit
------------
  The library is divided into five sub-modules: datastruct, mapping, sequencing
  testing, and util. After installing the libraries, they can be used either
  from the interactive python prompt or in scripts. To use, for example, the
  IntervalTree class which is in the intervalTree module in the datastruct
  package, type (either in the python interpretor or your script):
  ```python
    >>> from pyokit.datastruct.intervalTree import IntervalTree
  ```  

About Pyokit
------------
  What follows is a very brief description of the tools provided in Pyokit.

  ###Interface###
  Pyokit contains an extension to the built-in python command-line parser.
  Here's an example of using it to create a simple program with a command-line
  interface.
  ```python
    import os, sys
    from pyokit.interface.cli import CLI, Option

    def getUI():
      programName = os.path.basename(sys.argv[0])
      longDescription  =  "This script foos the biz baz and outputs a bar.\n" +\
                          "Usage: ./" + programName + " <biz.txt> <baz.txt>"
      shortDescription =  longDescription
      ui = CLI(programName, shortDescription, longDescription)
      ui.minArgs = 2
      ui.maxArgs = 2
      ui.addOption(Option(short = "v", long = "verbose",
                          description = "output status messages to stderr ",
                          required = False))
      ui.addOption(Option(short = "n", long = "num", argName = "fooNum",
                          description = "number of times to foo the biz baz",  
                          required = False, type = int))
      ui.addOption(Option(short = "o", long = "output_fn", argName = "filename",
                          description = "output bar to this file, else stdout",
                          required = False, type = str))
      ui.addOption(Option(short="h", long="help",
                          description="show this help message ", special=True))
      ui.parseCommandLine(sys.argv[1:])
      return ui

    def _main():
      ui = getUI()
      if ui.optionIsSet("help") :
        ui.usage()
        sys.exit()
      verbose = (ui.optionIsSet("verbose") == True)
      outfh = open(ui.getValue("output_fn")) \
              if ui.optionIsSet("output_fn") \
              else sys.stdout
      num = ui.getValue("num") if ui.optionIsSet("num") else 1
      biz = [l.strip() for l in open(ui.getArgument(0)) if l.strip() != ""]
      baz = [l.strip() for l in open(ui.getArgument(1)) if l.strip() != ""]
      for i in range(0, num) :
        if (verbose) : sys.stderr.write("fooing the biz baz\n")
        outfh.write("bar " + str(i+1) + ": " +\
                    ",".join(biz) + " -- " + ",".join(baz) + "\n")

    if __name__ == "__main__":
      _main()
  ```

  ###Iterators###
  Pyokit contains iterators for BED files, Wig files, MAF files, fasta files and
  fastq files.

  ***BED iterators***

  ***WIG iterators***

  ***MAF iterators***

  ***fasta iterators***

  ***fastq iterators***

  ###Bulk loading###
  You can easily load any of the above file types in full using the following
  syntax.

  In addition, Pyokit has some convenience functions for loading genomic regions
  into IntervalTree.

  ###Data structures###

Contacts and bug reports
------------------------
Philip J. Uren
philip.uren@gmail.com

If you found a bug or mistake in this project, I would like to know about it.
Before you send me the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been
   fixed.
2. Check that your input is in the correct format and you have selected the
   correct options.
3. Please reduce your input to the smallest possible size that still produces
   the bug; I will need your input data to reproduce the problem, and the
   smaller you can make it, the easier it will be for both of us.

License information
-------------------

  The following applys to this software package and all subparts therein  

Pyokit Copyright (C) 2014 Philip J. Uren

This library is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

This library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License along
with this library; if not, write to the Free Software Foundation, Inc., 51
Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
