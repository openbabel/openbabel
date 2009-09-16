@echo off

set OPTIONS=
set OPTIONS=-DMINIMAL_BUILD=ON

if NOT EXIST build mkdir build
cd build
"C:\Program Files\CMake 2.6\bin\cmake.exe" -DCMAKE_INSTALL_PREFIX=..\install -G "Visual Studio 9 2008" -DLIBXML2_LIBRARIES=..\..\windows-vc2005\libxml2.lib -DLIBXML2_INCLUDE_DIR=. -DZLIB_LIBRARY=..\..\windows-vc2005\zlib1.lib -DZLIB_INCLUDE_DIR=. %OPTIONS% ..\..
cd ..
