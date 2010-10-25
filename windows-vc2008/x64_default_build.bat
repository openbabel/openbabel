@rem -DMINIMAL_BUILD=OFF -DBUILD_GUI=ON -DENABLE_TESTS=ON

@if NOT EXIST build_x64 mkdir build_x64
@cd build_x64
"C:\Program Files (x86)\CMake 2.6\bin\cmake.exe" -DCMAKE_INSTALL_PREFIX=..\install -G "Visual Studio 9 2008 Win64" -DLIBXML2_LIBRARIES="%CD%\..\libs\x64\libxml2.lib" -DMINIMAL_BUILD=OFF -DWITH_INCHI=ON -DBUILD_GUI=OFF -DLIBXML2_INCLUDE_DIR=. -DZLIB_LIBRARY="%CD%\..\libs\x64\zlib1.lib" -DZLIB_INCLUDE_DIR=. -DINCHI_LIBRARY="%CD%\..\libs\x64\libstdinchi.lib" -DENABLE_TESTS=ON %1 %2 %3 %4 %5 %6 %7 %8 ..\..

@cd ..
