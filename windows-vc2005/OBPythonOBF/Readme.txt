Notes on compiling a Python interface to OpenBabel
using Visual Studio Express in Windows.

Python is assumed to be installed. 
Add the following environment variables
(Settings/Control Panel/System/Advanced/Environment Variables)
with either User or System scope. The values are typical examples.

PYTHON_DIR   C:\Python25

Install SWIG. Actually all that is needed is swig.exe which is
downloaded together with the source files and documentation.
Add its folder name to one of the PATH environment variables.

The project OBPythonOBF does the following:
 - It runs SWIG on the interface file scripts/openbable-python.i
    producing openbabel-python_wrap.cpp and openbabel.py
 - openbabel-python_wrap.cpp is compiled and linked with OBConv.lib, OBDLL.lib,
    and some other OB libs, to produce the DLL _openbabel.pyd
 - It runs obpython.bat which creates a Python distribution in dist
 
There is only a release configuration.
The SWIG step is not currently automatically rerun when any of the
OpenBabel files changes. If they have, you need to Rebuild the
OBPythonOBF project, which will also rebuild several OpenBabel DLLs first.

A useful article on building C++ extensions for Python is at
http://www.geocities.com/foetsch/python/extending_python.htm

Original notes:
 To use OpenBabel with Python either 
  run Python with the windows-vc2005 folder as the current folder, or
  add the windows-vc2005 folder to the PYTHONPATH environment variable.

New notes:
 To create a binary installer for the Python extension, do as follows...
 In the OBPythonOBF folder, run
     python setup.py bdist_wininst
                     --install-script=openbabel_postinstall.py
                     --bitmap=logo.bmp
 To install, run the created "openbabel-python-1.0.win32.exe" in the dist subfolder.
 You can test by copying testpybel.py to a folder that doesn't contain "openbabel.py*",
 and running it (this ensures that it uses the globally installed Open Babel;
 you should also ensure that PYTHONPATH does not contain a folder with "openbabel.py*").
 You can uninstall using Add/Remove Programs. 
 There is a log of the installation in $PYTHONDIR/openbabel-python-wininst.log.
