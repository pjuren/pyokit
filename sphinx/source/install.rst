=================
Installing Pyokit
=================

------------
Dependencies
------------
Pyokit depends on rpy2 for some of its functionality. You don't strictly
need it, Pyokit can be installed without it, but if you don't have it some
functionality will not be supported.

--------------------
Installing from PyPI
--------------------
Pyokit is available from the Python package index. The easiest
way to install it is to use pip. Note that You may need to prefix certain
commands below with ``sudo`` if you're installing to the global python
installation.

1. Get pip (if you don't already have it):

   .. code-block:: bash

     wget https://bootstrap.pypa.io/get-pip.py
     python get-pip.py

   This will install pip for you. if ``python`` doesn't point to the version of
   python you want to use, then replace it with the full path to the one you
   prefer (you'll probably have to do that if you don't have admin. access on
   the machine you're using).

2. Now you can install pyokit:

   .. code-block:: bash

     pip install pyokit

   ... and you're done.

----------------------
Installing from source
----------------------
The bleeding-edge version can always be cloned from the GitHub repository
(www.github.com/pjuren/pyokit), where tagged releases can also be downloaded as
.tar.gz files. If you take this route, you will still need pip (see above if you
don't have it). After unpacking the distribution (or cloning the repository),
``cd`` into the newly created directory and then

1. Build the distribution package:

   .. code-block:: bash

     python setup.py sdist

2. Install the distribution:

   .. code-block:: bash

     pip install --no-index dist/pyokit-a.b.c.tar.gz

   where a.b.c is the version you downloaded.
