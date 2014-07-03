# The following applies to this software package and all subparts therein
#
# Pyokit Copyright (C) 2014 Philip J. Uren
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this library; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA


################################################################################
##                              INSTALLATION                                  ##
################################################################################

uninstall :
	pip uninstall pyokit

install :
	rm -rf dist
	python setup.py sdist
	pip install --no-index dist/pyokit-?.?.?.tar.gz

################################################################################
##                              DOCUMENTATION                                 ##
################################################################################

docs :
	$(MAKE) -C sphinx clean
	$(MAKE) -C sphinx html
	rm -rf Docs/doctrees
	mv Docs/html/* Docs
	mv Docs/html/.buildinfo Docs
	rmdir Docs/html
.PHONY : docs

################################################################################
##                               UNIT TESTS                                   ##
################################################################################

test :
	python src/test.py

################################################################################
##                              HOUSEKEEPING                                  ##
################################################################################

clean :
	rm -f `find . -name "*.pyc"`

################################################################################
##                      VERSION AND RELEASE MANAGEMENT                        ##
################################################################################

releasePatch :
	bumpversion --verbose patch --tag --tag-name pyokit_{new_version}
.PHONY : releasePatch

releaseMinor :
	bumpversion --verbose minor --tag --tag-name pyokit_{new_version}
.PHONY : releaseMinor

releaseMajor :
	bumpversion --verbose major --tag --tag-name pyokit_{new_version}
.PHONY : releaseMajor

publishPyPITest :
	python setup.py sdist upload -r pypitest
.PHONY : publishPyPITest

publishPyPI :
	python setup.py sdist upload -r pypi
.PHONY : publishPyPI
