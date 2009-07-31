project(openbabel)
set(CMAKE_MODULE_PATH ${openbabel_SOURCE_DIR}/cmake/modules)
cmake_minimum_required(VERSION 2.4.8)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

ENABLE_TESTING()

include (CheckCXXCompilerFlag)
include (MacroEnsureVersion)

# Version numbering - should be bumped for each release
set(BABEL_MAJ_VER   2)
set(BABEL_MIN_VER   2)
set(BABEL_PATCH_VER 3)
if(WIN32)
  set(BABEL_LIBRARY   openbabel-2)
else(WIN32)
  set(BABEL_LIBRARY   openbabel)
endif(WIN32)
set(BABEL_VERSION  "${BABEL_MAJ_VER}.${BABEL_MIN_VER}.${BABEL_PATCH_VER}")
set(BABEL_DATADIR  "${CMAKE_INSTALL_PREFIX}/share/openbabel")

option(ENABLE_VERSIONED_FORMATS
  "Enable versioning of the format plugin directory" ON)

# Set up our path handling, inspired by the path handling used in KDE
set(BIN_INSTALL_DIR      "bin"
    CACHE PATH "Install dir for binaries")
set(LIB_SUFFIX "" CACHE STRING "Suffix of the directory name, e.g. 64 for lib64")
if(UNIX)
  set(LIB_INSTALL_DIR    "lib${LIB_SUFFIX}"
      CACHE PATH "Install dir for libraries")
  if(ENABLE_VERSIONED_FORMATS)
    set(OB_PLUGIN_INSTALL_DIR "${LIB_INSTALL_DIR}/openbabel/${BABEL_VERSION}"
      CACHE PATH "Install dir for plugins")
  else(ENABLE_VERSIONED_FORMATS)
    set(OB_PLUGIN_INSTALL_DIR "${LIB_INSTALL_DIR}/openbabel"
      CACHE PATH "Install dir for plugins")
  endif(ENABLE_VERSIONED_FORMATS)
else(UNIX) # Windows - bin dir = lib dir to load libraries
  set(LIB_INSTALL_DIR    "${BIN_INSTALL_DIR}"
      CACHE PATH "Install dir for libraries")
  if(ENABLE_VERSIONED_FORMATS)
    set(OB_PLUGIN_INSTALL_DIR "lib${LIB_SUFFIX}/openbabel/${BABEL_VERSION}"
      CACHE PATH "Install dir for plugins")
  else(ENABLE_VERSIONED_FORMATS)
    set(OB_PLUGIN_INSTALL_DIR "lib${LIB_SUFFIX}/openbabel"
      CACHE PATH "Install dir for plugins")
  endif(ENABLE_VERSIONED_FORMATS)
endif(UNIX)

option(WITH_INCHI   "build inchi lib" ON)
option(BUILD_SHARED "enable shared build support" ON)
if(BUILD_SHARED)
  set(BUILD_TYPE SHARED)
  set(PLUGIN_TYPE MODULE)
  add_definitions(-DUSING_DYNAMIC_LIBS)
else(BUILD_SHARED)
  set(BUILD_TYPE STATIC)
  set(PLUGIN_TYPE STATIC)
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
check_include_file(strings.h   HAVE_STRINGS_H)
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

   set(OB_MODULE_PATH "${CMAKE_INSTALL_PREFIX}/${OB_PLUGIN_INSTALL_DIR}")
   add_definitions(-DOB_MODULE_PATH="\\"${OB_MODULE_PATH}\\"")
endif(NOT WIN32)

# Add some visibility support when using GCC
if (CMAKE_COMPILER_IS_GNUCXX)
  if (CMAKE_SYSTEM_NAME MATCHES Linux)
    set (CMAKE_SHARED_LINKER_FLAGS "-Wl,--fatal-warnings -Wl,--no-undefined -lc ${CMAKE_SHARED_LINKER_FLAGS}")
    set (CMAKE_SHARED_LINKER_FLAGS "-Wl,--enable-new-dtags ${CMAKE_SHARED_LINKER_FLAGS}")
    set (CMAKE_MODULE_LINKER_FLAGS "-Wl,--fatal-warnings -Wl,--no-undefined -lc ${CMAKE_SHARED_LINKER_FLAGS}")
    set (CMAKE_MODULE_LINKER_FLAGS "-Wl,--enable-new-dtags ${CMAKE_SHARED_LINKER_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "-Wl,--fatal-warnings -Wl,--no-undefined -lc ${CMAKE_SHARED_LINKER_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-new-dtags ${CMAKE_EXE_LINKER_FLAGS}")
  endif (CMAKE_SYSTEM_NAME MATCHES Linux)
  # Now check if we can use visibility to selectively export symbols
  # Get the GCC version - from KDE4 cmake files
  exec_program(${CMAKE_C_COMPILER} ARGS --version OUTPUT_VARIABLE _gcc_version_info)
  string (REGEX MATCH "[345]\\.[0-9]\\.[0-9]" _gcc_version "${_gcc_version_info}")
  # gcc on mac just reports: "gcc (GCC) 3.3 20030304 ..." without the
  # patch level, handle this here:
  if (NOT _gcc_version)
    string (REGEX REPLACE ".*\\(GCC\\).* ([34]\\.[0-9]) .*" "\\1.0" _gcc_version "${_gcc_version_info}")
  endif (NOT _gcc_version)
  check_cxx_compiler_flag(-fvisibility=hidden HAVE_GCC_VISIBILITY)
  set(HAVE_GCC_VISIBILITY ${HAVE_GCC_VISIBILITY} CACHE BOOL "GCC support for hidden visibility")
  macro_ensure_version("4.1.0" "${_gcc_version}" GCC_IS_NEWER_THAN_4_1)
  macro_ensure_version("4.2.0" "${_gcc_version}" GCC_IS_NEWER_THAN_4_2)
  set(_GCC_COMPILED_WITH_BAD_ALLOCATOR FALSE)
  if (GCC_IS_NEWER_THAN_4_1)
    exec_program(${CMAKE_C_COMPILER} ARGS -v OUTPUT_VARIABLE _gcc_alloc_info)
    string(REGEX MATCH "(--enable-libstdcxx-allocator=mt)" _GCC_COMPILED_WITH_BAD_ALLOCATOR "${_gcc_alloc_info}")
  endif (GCC_IS_NEWER_THAN_4_1)
  if (HAVE_GCC_VISIBILITY AND GCC_IS_NEWER_THAN_4_1 AND NOT _GCC_COMPILED_WITH_BAD_ALLOCATOR)
    # We have all the parts necessary - use visibility support
    add_definitions(-DHAVE_GCC_VISIBILITY)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
    if (GCC_IS_NEWER_THAN_4_2)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility-inlines-hidden")
    endif (GCC_IS_NEWER_THAN_4_2)
  endif (HAVE_GCC_VISIBILITY AND GCC_IS_NEWER_THAN_4_1 AND NOT _GCC_COMPILED_WITH_BAD_ALLOCATOR)
endif (CMAKE_COMPILER_IS_GNUCXX)

# Some rpath handling for Linux and Mac
if(UNIX)
  if(APPLE)
    set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR})
  else(APPLE)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif(APPLE)
endif(UNIX)

# this is adapted from KDE's FindKDE4Internal.cmake : default the build type to
# "release with debug info".
#
# We will define two other build types: Debug and Release.
# These names are case-insensitive i.e. you can do -DCMAKE_BUILD_TYPE=debug
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif (NOT CMAKE_BUILD_TYPE)

#define various build types
if (CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
  set (CMAKE_CXX_FLAGS_DEBUG   "-g3 -fno-inline")
endif (CMAKE_COMPILER_IS_GNUCXX)

file(GLOB_RECURSE headers include/openbabel/*.h)
configure_file(${openbabel_SOURCE_DIR}/src/config.h.cmake ${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h)
install(FILES ${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h DESTINATION include/openbabel-2.0/openbabel)
configure_file(${openbabel_SOURCE_DIR}/openbabel-2.0.pc.cmake ${openbabel_BINARY_DIR}/openbabel-2.0.pc @ONLY)
install(FILES ${openbabel_BINARY_DIR}/openbabel-2.0.pc DESTINATION ${LIB_INSTALL_DIR}/pkgconfig)

include_directories(${openbabel_BINARY_DIR}/include
                    ${openbabel_SOURCE_DIR}/data
                    ${openbabel_SOURCE_DIR}/include
)

add_subdirectory(include)
add_subdirectory(data)
add_subdirectory(doc)
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(tools)

# Should Python bindings be built?
option(ENABLE_PYTHON "Enable Python bindings" OFF)
if(ENABLE_PYTHON)
  add_subdirectory(scripts)
endif(ENABLE_PYTHON)

