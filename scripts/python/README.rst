Open Babel Python Bindings
==========================

This is a Python interface to the Open Babel chemistry library. 

Open Babel is a chemical toolbox designed to speak the many languages
of chemical data. It's an open, collaborative project allowing anyone
to search, convert, analyze, or store data from molecular modeling,
chemistry, solid-state materials, biochemistry, or related areas.

For installation instructions, tutorials and examples, please visit the
`Open Babel website <http://openbabel.org/wiki/Python>`_.

This package provides two Python modules that can be used to access the
functionality of Open Babel toolkit:

-  The `openbabel <http://openbabel.org/docs/dev/UseTheLibrary/PythonDoc.html>`_
   module: A wrapper that is automatically generated using the SWIG package
   and provides access to almost all of the Open Babel interfaces via Python,
   including the base classes OBMol, OBAtom, OBBond, and OBResidue, as well 
   as the conversion framework OBConversion.

-  The `pybel <http://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html>`_
   module: A lightweight wrapper around the classes and methods in the 
   openbabel module. Pybel provides more convenient and Pythonic ways to 
   access the Open Babel toolkit.

Dependencies
------------

-  Python 2.4 or a more recent version.
-  Open Babel 2.3.0 or a more recent version.

Installation
------------

**Option 1**: Use `pip <http://www.pip-installer.org/en/latest/>`_.

::

    pip install openbabel-python

**Option 2**: Download the latest release and install yourself.

::

    tar -xzvf openbabel-python-1.8.tar.gz
    cd openbabel-1.8
    python setup.py install
    
**Option 3**: While building Open Babel itself.

::

    cd openbabel
    mkdir ob-build
    cd ob-build
    cmake -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON ..
    make install

Copyright and Licence
---------------------

-  Copyright (C) 2005-2007 Geoffrey R. Hutchison babel@geoffhutchison.net
-  Some portions Copyright (C) 2006-2010 Noel O'Boyle

This Python module is part of the `Open Babel project 
<http://openbabel.org/>`_.

Open Babel is distributed under the GNU General Public License (GPL).
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License. Full details
can be found in the file "COPYING" which should be included in your
distribution.
