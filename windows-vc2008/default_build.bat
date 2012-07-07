@rem The cmake/bin directory should be in PATH
@rem To build the GUI, the top level directory of wxWidgets should be in the
@rem environment variable WXWIN and -DBUILD_GUI=ON specified 
@if NOT EXIST build mkdir build
@cd build
cmake.exe -DCMAKE_INSTALL_PREFIX=..\install -G "Visual Studio 10" -DLIBXML2_LIBRARIES="%CD%\..\libs\i386\libxml2.lib" -DMINIMAL_BUILD=OFF -DCAIRO_INCLUDE_DIRS="%CD%\..\include\cairo" -DCAIRO_LIBRARIES="%CD%\..\libs\i386\cairo.lib" -DLIBXML2_INCLUDE_DIR=. -DZLIB_LIBRARY="%CD%\..\libs\i386\zlib1.lib" -DZLIB_INCLUDE_DIR=. -DINCHI_LIBRARY="%CD%\..\libs\i386\libinchi.lib" -DINCHI_INCLUDE_DIR=. -DXDR_LIBRARY="%CD%\..\libs\i386\xdr.lib" -DRUN_SWIG=ON -DJAVA_BINDINGS=ON -DCSHARP_BINDINGS=ON -DCSHARP_EXECUTABLE=C:\Windows\Microsoft.NET\Framework\v3.5\csc.exe -DENABLE_TESTS=OFF -DBUILD_GUI=ON -DwxWidgets_ROOT_DIR="%WXWIN%" -DwxWidgets_LIB_DIR="%WXWIN%/lib/vc_lib" -DwxWidgets_CONFIGURATION=msw %1 %2 %3 %4 %5 %6 %7 %8 ..\..
@cd ..
