@rem -DMINIMAL_BUILD=OFF -DBUILD_GUI=ON -DENABLE_TESTS=ON

@rem The cmake/bin directory should be in PATH
@rem The top level directory of Eigen should be in the the environment variable
@rem EIGEN2_INCLUDE_DIR. If it is not found then qeq.cpp and qtpie.cpp are not built.
@if NOT EXIST build mkdir build
@cd build
cmake.exe -DCMAKE_INSTALL_PREFIX=..\install -G "Visual Studio 9 2008" -DLIBXML2_LIBRARIES="%CD%\..\libs\i386\libxml2.lib" -DMINIMAL_BUILD=OFF -DBUILD_GUI=ON -DLIBXML2_INCLUDE_DIR=. -DZLIB_LIBRARY="%CD%\..\libs\i386\zlib1.lib" -DZLIB_INCLUDE_DIR=. -DINCHI_LIBRARY="%CD%\..\libs\i386\libinchi.lib" -DXDR_LIBRARY="%CD%\..\libs\i386\xdr.lib" -DEIGEN2_INCLUDE_DIR="%EIGEN2_INCLUDE_DIR%" -DENABLE_TESTS=OFF %1 %2 %3 %4 %5 %6 %7 %8 ..\..
@cd ..
