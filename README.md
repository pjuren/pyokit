Pyokit: Python Bio Toolkit
==========================
[![Build Status](https://travis-ci.org/pjuren/pyokit.svg?branch=master)](https://travis-ci.org/pjuren/pyokit)
[![PyPI version](https://badge.fury.io/py/pyokit.svg)](http://badge.fury.io/py/pyokit)
[![Code Health](https://landscape.io/github/pjuren/pyokit/master/landscape.svg?style=flat)](https://landscape.io/github/pjuren/pyokit/master)
[![Coverage Status](https://coveralls.io/repos/pjuren/pyokit/badge.svg)](https://coveralls.io/r/pjuren/pyokit)

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
way to install it is to use pip. Note that you may need to prefix certain
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
Everything you need to know should be in the Pyokit manual. You can find
it here: http://pjuren.github.io/pyokit

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
