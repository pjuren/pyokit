#  Copyright (C) 2011
#  University of Southern California,
#  Philip J. Uren,
#  
#  Authors: Philip J. Uren
#  
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  


# Change this to point to the python install where you
# want to install the library
PYTHON_PATH=python

# build the python egg
build :
	cd src;	$(PYTHON_PATH) setup.py bdist_egg

# install the python egg
install : 
	cd src;	$(PYTHON_PATH) setup.py install

clean :
	cd src; rm -rf smithlab_py.egg-info build dist

test :
	cd src; $(PYTHON_PATH) test.py

all : build
