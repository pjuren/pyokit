
# build the python egg
build :
	cd src;	python setup.py bdist_egg

# install the python egg
install : 
	cd src;	python setup.py install

clean :
	cd src; rm -rf smithlab_py.egg-info build dist

all : build
