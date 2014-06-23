import sys, os

# in case the user doesn't have setup tools, or the right version
import ez_setup
ez_setup.use_setuptools()

# should be able to safely do this now.
from setuptools import setup, find_packages

setup(name='pyokit',
      version='0.1',
      packages = find_packages('src'),  # include all packages under src
			package_dir = {'':'src'},   # ell distutils packages are under src
			description = 'A library of python functions and classes to assist in processing high-throughput biological data',
			author = 'Philip J. Uren',
			author_email = 'philip.uren@gmail.com',
			url = 'https://github.com/pjuren/pyokit', 
			download_url = 'https://github.com/pjuren/pyokit/tarball/pyokit_0.1.1', 
			keywords = [],
			classifiers = [],
)

