# change this to point to the python install where you
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
