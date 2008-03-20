project(openbabel)
#set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules )
cmake_minimum_required(VERSION 2.4.5)

set(BABEL_MAJ_VER   2)
set(BABEL_MIN_VER   2)
set(BABEL_PATCH_VER 0b4)
if(WIN32)
  set(BABEL_LIBRARY   openbabel-2)
else(WIN32)
  set(BABEL_LIBRARY   openbabel)
endif(WIN32)
set(BABEL_VERSION  "${BABEL_MAJ_VER}.${BABEL_MIN_VER}.${BABEL_PATCH_VER}")
set(BABEL_DATADIR  "${CMAKE_INSTALL_PREFIX}/share/openbabel")
set(LIB            "lib")   ## FIXME: add lib suffix!

option(WITH_INCHI   "build inchi lib" ON)
option(BUILD_SHARED "enable static build support" ON)
if(BUILD_SHARED)
  set(BUILD_TYPE SHARED)
  add_definitions(-DUSING_DYNAMIC_LIBS)
else(BUILD_SHARED)
  set(BUILD_TYPE STATIC)
endif(BUILD_SHARED)

if(WIN32)
  # FIXME: not used on windows... - inconsistency
  set(MODULE_EXTENSION ".obf")
else(WIN32)
  set(MODULE_EXTENSION ".so")
endif(WIN32)

# configure checks
find_package(LibXml2)
find_package(Boost)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
  add_definitions(-DUSE_BOOST)
endif(Boost_FOUND)
find_package(ZLIB)
if(ZLIB_FOUND)
  add_definitions(-DHAVE_LIBZ)
  include_directories(${ZLIB_INCLUDE_DIR})
endif(ZLIB_FOUND)

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckSymbolExists)
include(CheckTypeSize)
include(CheckCSourceCompiles)

check_include_file(conio.h      HAVE_CONIO_H)
check_include_file(sys/time.h   HAVE_SYS_TIME_H)
check_include_file(time.h       HAVE_TIME_H)
check_include_file_cxx(sstream  HAVE_SSTREAM)

check_symbol_exists(rint          "math.h"     HAVE_RINT)
check_symbol_exists(snprintf      "stdio.h"    HAVE_SNPRINTF)
check_symbol_exists(sranddev      "stdlib.h"   HAVE_SRANDDEV)
check_symbol_exists(strcasecmp    "string.h"   HAVE_STRCASECMP)
check_symbol_exists(strncasecmp   "string.h"   HAVE_STRNCASECMP)

set(CMAKE_EXTRA_INCLUDE_FILES time.h)
check_type_size(clock_t CLOCK_T)

if(NOT WIN32)
  SET(CMAKE_REQUIRED_FLAGS -Werror)
  check_c_source_compiles("
   #include <sys/types.h>
   #include <dirent.h>
    int main(){
     extern int matchFiles (struct dirent *entry_p);
     struct dirent **entries_pp;
     int count = scandir (\"./\", &entries_pp, matchFiles, 0);
     return 0;
    }
   " SCANDIR_NEEDS_CONST)

   set(OB_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/${LIB}/openbabel)
endif(NOT WIN32)

file(GLOB_RECURSE headers include/openbabel/*.h)
configure_file(${CMAKE_SOURCE_DIR}/src/config.h.cmake ${CMAKE_BINARY_DIR}/include/openbabel/babelconfig.h)
install(FILES ${CMAKE_BINARY_DIR}/include/openbabel/babelconfig.h DESTINATION include/openbabel-2.0/openbabel)
configure_file(${CMAKE_SOURCE_DIR}/openbabel-2.0.pc.cmake ${CMAKE_BINARY_DIR}/openbabel-2.0.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/openbabel-2.0.pc DESTINATION ${LIB}/pkgconfig)

include_directories(${CMAKE_SOURCE_DIR}/include
                    ${CMAKE_SOURCE_DIR}/data
                    ${CMAKE_BINARY_DIR}/include
)

add_subdirectory(include)
add_subdirectory(data)
add_subdirectory(doc)
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(tools)

