@rem The cmake/bin directory should be in PATH
@rem The top level directory of Eigen should be in the the environment variable
@rem EIGEN2_INCLUDE_DIR.
@if NOT EXIST build_x64 mkdir build_x64
@cd build_x64
cmake.exe -DCMAKE_INSTALL_PREFIX=..\install -G "Visual Studio 10 Win64" -DLIBXML2_LIBRARIES="%CD%\..\libs\x64\libxml2.lib" -DWITH_INCHI=ON -DHAVE_RPC_XDR_H=FALSE -DMINIMAL_BUILD=OFF -DBUILD_GUI=OFF -DLIBXML2_INCLUDE_DIR=. -DZLIB_LIBRARY="%CD%\..\libs\x64\zlib1.lib" -DZLIB_INCLUDE_DIR=. -DINCHI_LIBRARY="%CD%\..\libs\x64\libinchi.lib" -DINCHI_INCLUDE_DIR=. -DEIGEN2_INCLUDE_DIR="%EIGEN2_INCLUDE_DIR%" -DRUN_SWIG=ON -DJAVA_BINDINGS=ON -DENABLE_TESTS=OFF %1 %2 %3 %4 %5 %6 %7 %8 ..\..
@cd ..
