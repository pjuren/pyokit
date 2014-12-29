=======================
Developer documentation
=======================


--------------------
Testing your changes
--------------------
All of the test scripts (and indeed anything else you write that uses Pyokit)
are using the INSTALLED version of the library, not the source. So if you
make changes, you'll need to rebuild and reinstall the library to test them
out.

The Python package manager will not install Pyokit again if the version number
has not changed -- which it often won't while you're making changes. To get
around this, you just need to uninstall the library first, the install it again.
So, in practice, this:

.. code-block:: bash

  make uninstall; make install

-------------------------
Maintaining documentation
-------------------------
Pyokit uses Sphinx for documentation. All of these docs are constructed
automatically from the Pyokit source code and the documentation source
(which can be found in the sphinx/source directory). The documentation is
hosted on Github pages. Travis CI is used to automatically build it and
push the new copy every time there is a push to the Pyokit repository on
GitHub.

To update the documentation, you should edit the sphinx source .rst files.
The Makefile has a target that will rebuild the (local) documentation for you:

.. code-block:: bash

  make docs

Sphinx will tell you where the output documentation was placed (Docs/html at
present, relative to the project root). These are just HTML files and you
can open the Index.html in any web browser to preview the documentation.

.. warning:: It is the INSTALLED version of Pyokit source code that is used
             when building these, but the SOURCE .rst files. Practically, that
             means you can get away with just rebuilding the documentation if
             you are only editing the Sphinx .rst source files, but changes
             to the Pyokit source will not register until you rebuild and
             install the library itself locally.
