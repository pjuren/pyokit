
clean :
	rm -f `find . -name "*.pyc"`

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

releasePatch :
	bumpversion --verbose patch --tag --tag-name pyokit_{new_version}
.PHONY : releasePatch

releaseMinor :
	bumpversion --verbose minor --tag --tag-name pyokit_{new_version}
.PHONY : releaseMinor

releaseMajor :
	bumpversion --verbose major --tag --tag-name pyokit_{new_version}
.PHONY : releaseMajor

