OpenBabel Version 2.0.2 Build files for Windows

The project files to build the OpenBabel command line program are provided for
for two versions of Visual Studio version 6, and version 8. For the latter,
Visual Studio 2005 Express (available free until Nov 2006) is sufficient.

Open OpenBabel.dsw (VC6) or OpenBabel.sln (VC8) in Visual Studio and build
the OBabel project in "Release" build. They can be also built in "Debug" build,
but for VC6 the very large number of warnings about variable name length can be
ignored.

Any of the cpp files implementing the formats can be excluded from the the
project without any other modifications being necessary. To write your own format
see  the file exampleformat.cpp and 
http://openbabel.sourceforge.net/wiki/HowTo:Add_A_New_File_Format

OpenBabel was designed so that various modules could be separately compiled as
DLLs. The main API containing the chemistry, the conversion routines, 
any number of formats and user interfaces (command line or GUI) can be built
as DLLs or (for the interface) and exe. The project files for these builds are
not provided in this distribution because a few unexplained problems were
encountered. It is hoped to resolve these soon, and if you feel that this type
of build would be useful to you, please post a note to the development mailing list
(see http://openbabel.sourceforge.net/wiki/Main_Page) so that we can judge the
potential interest.

A VC8 project to build  the OpenBabelGUI with the OpenBabel version 2.0.1 source
files is also included. It is based on wxWidgets, which needs to be installed first.
This is will eventually be cross-platform and the Windows version here is currently
under development; building it may not be entirely straightforward. The compiled
version is usable out of the box.

     
