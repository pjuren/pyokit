
test : 
	python src/test.py

uninstall :
	pip uninstall pyokit

install :
	rm -rf dist
	python setup.py sdist
	pip install --no-index dist/pyokit-?.?.?.tar.gz

docs :
	$(MAKE) -C sphinx clean
	$(MAKE) -C sphinx html 
	rm -rf Docs/doctrees
	mv Docs/html/* Docs
	mv Docs/html/.buildinfo Docs
	rmdir Docs/html
.PHONY : docs
