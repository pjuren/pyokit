"""
Date of Creation: 2013.

Description:   Settings for setup tools to install this package

Copyright (C) 2010-2015
Philip J. Uren,

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
"""

from setuptools import setup, find_packages

setup(name='pyokit',
      version='0.2.0',
      packages=find_packages('src'),  # include all packages under src
      entry_points={'console_scripts':
                    ['pyokit=pyokit.scripts.pyokitMain:main']},
      package_dir={'': 'src'},   # ell distutils packages are under src
      install_requires=['mock>=1.0.0', 'enum34', 'subprocess32'],
      description='A library of python functions and classes to assist in '
                  'processing high-throughput biological data',
      author='Philip J. Uren',
      author_email='philip.uren@gmail.com',
      url='https://github.com/pjuren/pyokit',
      download_url='https://github.com/pjuren/pyokit/tarball/pyokit_0.2.0',
      license='LGPL2',
      keywords=[],
      classifiers=[])
