Pyokit: Python Bio Toolkit
==========================

Pyokit is a python library with a collection of tools to make the processing
and analysis of high-throughput biological datasets in python easier. 

Dependencies
------------
Pyokit depends on rpy2 for some of its functionality. You don't strictly
need it, Pyokit can be installed without it, but if you don't have it some 
functionality will not be supported. 

Configuring 
-----------
The package can be installed without any configuration using the instructions
in the section 'INSTALLING'. If you wish to install the package to a python
version/installation other than the default (i.e. the one invoked when 
typing 'python' on the command line) then you need to edit the Makefile and 
set the variable PYTHON_PATH to point to the desired python executable 
(in future this will be made easier...)

Installing
----------
This should be very easy. Just type the following at the command prompt,
assuming your current working directory is the root directory of the 
distribution:

  > make ; make install

Using Pyokit
------------
  The library is divided into five sub-modules: datastruct, mapping, sequencing
  testing, and util. After installing the libraries, they can be used either from
  the interactive python prompt or in scripts. To use, for example, the 
  IntervalTree class which is in the intervalTree module in the datastruct 
  package, type (either in the python interpretor or your script):
        
    >>> from pyokit.datastruct.intervalTree import IntervalTree

About Pyokit
------------
  What follows is a very breif description of the tools provided in Pyokit.
  For a more detailed explanation, please see the Pyokit manual. 

Contacts and bug reports
------------------------
Philip J. Uren philip.uren@gmail.com

If you found a bug or mistake in this project, I would like to know about it. 
Before you send me the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been fixed.
2. Check that your input is in the correct format and you have selected the correct options.
3. Please reduce your input to the smallest possible size that still produces the bug; 
   I will need your input data to reproduce the problem, and the smaller you can make it, 
   the easier it will be for both of us.

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

