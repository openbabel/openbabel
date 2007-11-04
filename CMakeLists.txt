project(openabel)
set(BABEL_MAJ_VER   2)
set(BABEL_MIN_VER   1)
set(BABEL_PATCH_VER 1)
set(BABEL_VERSION  "${BABEL_MAJ_VER}.${BABEL_MIN_VER}.${BABEL_PATCH_VER}")
set(BABEL_DATADIR  "{CMAKE_INSTALL_PREFIX}/share/openbabel")

option(WITH_INCHI   "build inchi lib" ON)
option(BUILD_SHARED "enable static build support" ON)
if(BUILD_SHARED)
  set(BUILD_TYPE SHARED)
else(BUILD_SHARED)
  set(BUILD_TYPE STATIC)
endif(BUILD_SHARED)

if(WIN32)
  set(MODULE_EXTENSION ".dll")
else(WIN32)
  set(MODULE_EXTENSION ".so")
endif(WIN32)

# configure checks
find_package(LibXml2)
find_package(ZLIB)
if(ZLIB_FOUND)
  add_definitions(-DHAVE_LIBZ)
  include_directories(${ZLIB_INCLUDE_DIR})
endif(ZLIB_FOUND)

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckSymbolExists)
include(CheckTypeSize)

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

configure_file(${CMAKE_SOURCE_DIR}/src/config.h.cmake ${CMAKE_BINARY_DIR}/include/openbabel/babelconfig.h)
install(FILES ${CMAKE_BINARY_DIR}/include/openbabel/babelconfig.h DESTINATION include/openbabel-2.0/openbabel)

include_directories(${CMAKE_SOURCE_DIR}/include
                    ${CMAKE_SOURCE_DIR}/data
                    ${CMAKE_BINARY_DIR}/include
)

add_subdirectory(include)
add_subdirectory(data)
add_subdirectory(src)

