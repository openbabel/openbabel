xdr.lib and xdr-0.dll are the import library and shared build of bsd-xdr-1.0.0 available from http://code.google.com/p/bsd-xdr/. These were built from source by starting the MSVC2008 shell, starting bash (provided by Cygwin), and running the command "make -f Makefile.msvc80".

The include files in %ROOT%/windows-vc2008/include are also from this code base.

Cairo and its dependencies (freetype, expat, fontconfig, libpng) were downloaded from 
http://www.gtk.org/download-windows.html. The Cairo include files and cairo.lib were taken
from Cairo Dev from that same webpage.
