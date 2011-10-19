================================================================================
                     SMITHLAB PYTHON COMMON LIBRARIES 
================================================================================

LICENSE
-------
  The following applys to this software package and all subparts therein  

  Copyright (C) 2011
  University of Southern California,
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
  

CONFIGURING
-----------
  The package can be installed without any configuration using the instructions
  in the section 'INSTALLING'. If you wish to install the package to a python
  version/installation other than the default (i.e. the one invoked when 
  typing 'python' on the command line) then you need to edit the Makefile and
  set the variable PYTHON_PATH to point to the desired python executable

INSTALLING 
----------
  > make ; make install

USING 
-----
  The source is divided into five packages, datastruct, mapping, sequences, 
  testing, util. After installing the libraries, they can be used either from
  the interactive python prompt or in scripts. To use, for example, the 
  IntervalTree class which is in the intervalTree module in the datastruct 
  package, type:
    
    >>> from datastruct.intervalTree import IntervalTree
