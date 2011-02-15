# change this to point to the python install where you
# want to install the library
PYTHON_PATH=python

# build the python egg
build :
	cd src;	python setup.py bdist_egg

# install the python egg
install : 
	cd src;	python setup.py install

clean :
	cd src; rm -rf smithlab_py.egg-info build dist

test :
	cd src; for dir in `ls -d */`; do cd $$dir; echo `pwd`; for file in `ls *.py`; do echo $$file; python $$file; done; cd ..; done

all : build
