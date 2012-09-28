# Please ensure that any changes remain compliant with 2.4.8.
if(NOT EMBED_OPENBABEL)
  cmake_minimum_required(VERSION 2.4.8)
endif()

project(openbabel)
set(CMAKE_MODULE_PATH ${openbabel_SOURCE_DIR}/cmake/modules)

# Allow loose loop constructs, i.e. no matching in if/else/endif or loops.
# Note that this is true by default in CMake 2.6.0, but we currently only
# require CMake 2.4.8 - remove this when the CMake requirement is bumped.
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS 1)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0005 OLD) # add_definitions need updating to set to NEW
endif()

include (CheckCXXCompilerFlag)
include (MacroEnsureVersion)

# Version numbering - should be bumped for each release
# Note that for "beta" releases, we should start at x.90.0 -- we've
# had too many 1.100.1 releases. :-)
set(BABEL_MAJ_VER   2)
set(BABEL_MIN_VER   3)
set(BABEL_PATCH_VER 2)

# This should be phased out in preference for just using the target name
set(BABEL_LIBRARY openbabel)

# This would be better handled with external projects...
if(EMBED_OPENBABEL)
  set(BABEL_VERSION  "${BABEL_MAJ_VER}.${BABEL_MIN_VER}.${BABEL_PATCH_VER}" PARENT_SCOPE)
else()
  set(BABEL_VERSION  "${BABEL_MAJ_VER}.${BABEL_MIN_VER}.${BABEL_PATCH_VER}")
endif()
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
  else()
    set(OB_PLUGIN_INSTALL_DIR "${LIB_INSTALL_DIR}/openbabel"
      CACHE PATH "Install dir for plugins")
  endif()
else() # Windows - bin dir = lib dir to load libraries
  set(LIB_INSTALL_DIR "${BIN_INSTALL_DIR}"
      CACHE PATH "Install dir for libraries")
  if(MSVC)
    set(OB_PLUGIN_INSTALL_DIR "${BIN_INSTALL_DIR}"
      CACHE PATH "Install dir for plugins")
  else()
    if(ENABLE_VERSIONED_FORMATS)
      set(OB_PLUGIN_INSTALL_DIR "lib${LIB_SUFFIX}/openbabel/${BABEL_VERSION}"
        CACHE PATH "Install dir for plugins")
    else(ENABLE_VERSIONED_FORMATS)
      set(OB_PLUGIN_INSTALL_DIR "lib${LIB_SUFFIX}/openbabel"
        CACHE PATH "Install dir for plugins")
    endif(ENABLE_VERSIONED_FORMATS)
  endif()
endif()
set(OB_INCLUDE_DIRS "include/openbabel-2.0")
set(OB_EXPORTS_FILE "${openbabel_BINARY_DIR}/OpenBabel2_EXPORTS.cmake")
# Ensure a fresh file is made each time CMake is run
file(REMOVE "${OB_EXPORTS_FILE}")

# Place binaries and libraries in the root of the build directory
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    CACHE STRING "Binary build output directory")
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  if(UNIX)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib${LIB_SUFFIX}"
      CACHE STRING "Library build output directory")
  else() # This is preferable on WIndows for example - bin = lib directory.
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
      CACHE STRING "Library build output directory")
  endif()
endif()

mark_as_advanced(CMAKE_LIBRARY_OUTPUT_DIRECTORY
  CMAKE_RUNTIME_OUTPUT_DIRECTORY)

option(BUILD_SHARED "enable shared build support" ON)
option(BUILD_MIXED "enable linux relocatable binaries support" OFF)
option(WITH_STATIC_INCHI "enable static inchi and related plugins" OFF)
option(WITH_STATIC_LIBXML "enable static libxml and xml file format plugins" OFF)
# Static building
if(BUILD_SHARED)
  set(BUILD_TYPE SHARED)
  set(PLUGIN_TYPE MODULE)
  add_definitions(-DUSING_DYNAMIC_LIBS)
else()
  set(BUILD_TYPE STATIC)
  set(PLUGIN_TYPE STATIC)

  if(WITH_STATIC_INCHI)
    add_definitions(-DHAVE_STATIC_INCHI)
  endif()
  if(WITH_STATIC_LIBXML)
    add_definitions(-DHAVE_STATIC_LIBXML)
  endif()

  if(MSVC)
    # set cl flags for static compiling
    set(CMAKE_CXX_FLAGS_DEBUG "/MTd")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "/INCREMENTAL:NO /NODEFAULTLIB:MSVCRT")
    set(CMAKE_CXX_FLAGS_RELEASE	"/MT /O2 /Ob2 /D NDEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/INCREMENTAL:NO /NODEFAULTLIB:MSVCRT")
    # note: static libraries are specified when running cmake
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  endif()
endif()

# Should we attempt to use the system inchi library?
option(WITH_INCHI   "Build inchi support" ON)
option(OPENBABEL_USE_SYSTEM_INCHI
  "Use the system inchi library." OFF)
if(OPENBABEL_USE_SYSTEM_INCHI)
  find_package(Inchi REQUIRED)
  if (NOT INCHI_FOUND)
    message(FATAL_ERROR "Inchi library not found.")
  endif()
else()
  message(STATUS "Using included inchi library.")
  set(INCHI_LIBRARIES inchi)
endif()

if(WIN32)
  # FIXME: not used on windows... - inconsistency
  set(MODULE_EXTENSION ".obf")
else()
  set(MODULE_EXTENSION ".so")
endif()

# configure checks
find_package(LibXml2)
if(NOT LIBXML2_FOUND)
  message(WARNING "libxml2 not found - disabling CML support!")
endif()

find_package(ZLIB)
if(ZLIB_FOUND)
  add_definitions(-DHAVE_LIBZ)
  include_directories(${ZLIB_INCLUDE_DIR})
  # Longstanding unsolved problem with compression under Windows
  if(WIN32)
    add_definitions(-DDISABLE_WRITE_COMPRESSION)
  endif()
endif()

# wxWidgets instructions based on http://wiki.wxwidgets.org/CMake
#find_package(wxWidgets COMPONENTS base core REQUIRED)
find_package(wxWidgets COMPONENTS base core adv)
if(wxWidgets_FOUND)
  include(${wxWidgets_USE_FILE})
  add_definitions(-DHAVE_WXWIDGETS)
  include_directories(${wxWidgets_INCLUDE_DIRS})
endif()

if(MSVC)
  # Ensure that CharacterSet="0" in the project files
  add_definitions(-D_SBCS) # Single-Byte Character Set (requires CMake 2.8.8)
endif()

if(MSVC)
  option(OB_USE_PREBUILT_BINARIES
    "Should Windows Open Babel builds use prebuilt binaries?" ON)
	# dlhandler_win32.cpp assumes multibyte character set
	remove_definitions(-DUNICODE -D_UNICODE)
endif()

if(MSVC AND OB_USE_PREBUILT_BINARIES)
  include_directories(${openbabel_SOURCE_DIR}/windows-vc2008/include)
  # The following line is required for check_include_file(rpc/xdr.h ...) below
  set(CMAKE_REQUIRED_INCLUDES ${openbabel_SOURCE_DIR}/windows-vc2008/include)
  # Force OPENBABEL_USE_SYSTEM_INCHI to ON, as this should be using the binary
  set(OPENBABEL_USE_SYSTEM_INCHI ON CACHE BOOL
    "Forced to ON for prebuilt binaries" FORCE)
endif()

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckSymbolExists)
include(CheckLibraryExists)
include(CheckTypeSize)
include(CheckCSourceCompiles)

check_include_file(conio.h      HAVE_CONIO_H)
check_include_file(sys/time.h   HAVE_SYS_TIME_H)
check_include_file(time.h       HAVE_TIME_H)
check_include_file(strings.h    HAVE_STRINGS_H)
check_include_file(rpc/xdr.h    HAVE_RPC_XDR_H)
check_include_file(regex.h      HAVE_REGEX_H)
check_include_file_cxx(sstream  HAVE_SSTREAM)

check_symbol_exists(rint          "math.h"     HAVE_RINT)
check_symbol_exists(snprintf      "stdio.h"    HAVE_SNPRINTF)
check_symbol_exists(sranddev      "stdlib.h"   HAVE_SRANDDEV)
check_symbol_exists(strcasecmp    "string.h"   HAVE_STRCASECMP)
check_symbol_exists(strncasecmp   "string.h"   HAVE_STRNCASECMP)

# BSDs don't link against libdl, but rather libc
check_library_exists(dl dlopen "" HAVE_LIBDL)

set(CMAKE_EXTRA_INCLUDE_FILES time.h)
check_type_size(clock_t CLOCK_T)

# Get the GCC version - from KDE4 cmake files
if(CMAKE_COMPILER_IS_GNUCXX)
  exec_program(${CMAKE_C_COMPILER} ARGS --version OUTPUT_VARIABLE _gcc_version_info)
  string(REGEX MATCH "[345]\\.[0-9]\\.[0-9]" _gcc_version "${_gcc_version_info}")
  # gcc on mac just reports: "gcc (GCC) 3.3 20030304 ..." without the
  # patch level, handle this here:
  if (NOT _gcc_version)
    string (REGEX REPLACE ".*\\(GCC\\).* ([34]\\.[0-9]) .*" "\\1.0" _gcc_version "${_gcc_version_info}")
  endif ()
  macro_ensure_version("4.0.0" "${_gcc_version}" GCC_IS_NEWER_THAN_4_0)
  macro_ensure_version("4.1.0" "${_gcc_version}" GCC_IS_NEWER_THAN_4_1)
  macro_ensure_version("4.2.0" "${_gcc_version}" GCC_IS_NEWER_THAN_4_2)
endif()

if(UNIX)
  option(ENABLE_SYMBOL_VISIBILITY
    "Enable support for compiler symbol visibility. WARNING: Some issues reported"
    OFF)
endif()

if(NOT MSVC)
  SET(CMAKE_REQUIRED_FLAGS -Werror)
  check_c_source_compiles("
   #include <sys/types.h>
   #include <dirent.h>
    int main(){
     extern int matchFiles (struct dirent *entry_p);
     struct dirent **entries_pp;
     int count = scandir (\"./\", &entries_pp, matchFiles, 0);
     return count;
    }
   " SCANDIR_NEEDS_CONST)

  set(OB_MODULE_PATH "${CMAKE_INSTALL_PREFIX}/${OB_PLUGIN_INSTALL_DIR}")
  add_definitions(-DOB_MODULE_PATH="\\"${OB_MODULE_PATH}\\"")

  # Add some visibility support when using GCC
  # note: Altough MinGW g++ 4.4 passes this test, visibility can't be used
  if(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)
    if(CMAKE_SYSTEM_NAME MATCHES Linux)
      set (CMAKE_SHARED_LINKER_FLAGS "-Wl,--fatal-warnings -Wl,--no-undefined -lc ${CMAKE_SHARED_LINKER_FLAGS}")
      set (CMAKE_SHARED_LINKER_FLAGS "-Wl,--enable-new-dtags ${CMAKE_SHARED_LINKER_FLAGS}")
      set (CMAKE_MODULE_LINKER_FLAGS "-Wl,--fatal-warnings -Wl,--no-undefined -lc ${CMAKE_SHARED_LINKER_FLAGS}")
      set (CMAKE_MODULE_LINKER_FLAGS "-Wl,--enable-new-dtags ${CMAKE_SHARED_LINKER_FLAGS}")
      if(BUILD_SHARED)
        set (CMAKE_EXE_LINKER_FLAGS "-Wl,--fatal-warnings -Wl,--no-undefined -lc ${CMAKE_SHARED_LINKER_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-new-dtags ${CMAKE_EXE_LINKER_FLAGS}")
      endif()
    endif ()
    # Now check if we can use visibility to selectively export symbols
    check_cxx_compiler_flag(-fvisibility=hidden HAVE_GCC_VISIBILITY)
    set(HAVE_GCC_VISIBILITY ${HAVE_GCC_VISIBILITY} CACHE BOOL "GCC support for hidden visibility")
    set(_GCC_COMPILED_WITH_BAD_ALLOCATOR FALSE)
    if(GCC_IS_NEWER_THAN_4_1)
      exec_program(${CMAKE_C_COMPILER} ARGS -v OUTPUT_VARIABLE _gcc_alloc_info)
      string(REGEX MATCH "(--enable-libstdcxx-allocator=mt)" _GCC_COMPILED_WITH_BAD_ALLOCATOR "${_gcc_alloc_info}")
    endif()
    if(HAVE_GCC_VISIBILITY AND GCC_IS_NEWER_THAN_4_2 AND ENABLE_SYMBOL_VISIBILITY
       AND NOT _GCC_COMPILED_WITH_BAD_ALLOCATOR)
      # We have all the parts necessary - use visibility support
      add_definitions(-DHAVE_GCC_VISIBILITY)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility-inlines-hidden")
    endif()
  endif()
endif()

# To avoid warnings when linking the plugins and exes
if(MINGW)
  set (CMAKE_MODULE_LINKER_FLAGS "-Wl,--enable-auto-import ${CMAKE_MODULE_LINKER_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-auto-import ${CMAKE_EXE_LINKER_FLAGS}")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  if(NOT GCC_IS_NEWER_THAN_4_0) # GCC 3.x
    find_package(Boost)
    if(Boost_FOUND)
      include_directories(${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
      add_definitions(-DUSE_BOOST -DHAVE_SHARED_POINTER)
      set(SHARED_POINTER true)
    else()
      message(FATAL_ERROR "Boost is required to compile OpenBabel with GCC 3.x")
    endif()
  else()
    add_definitions(-DHAVE_SHARED_POINTER)
    set(SHARED_POINTER true)
  endif()
endif()

if (${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
   set(SHARED_POINTER true)
   add_definitions(-DHAVE_SHARED_POINTER)
endif()

# Intel C++ Compiler
# Note: Windows version uses different options!
if(${CMAKE_CXX_COMPILER_ID} MATCHES Intel AND UNIX)
  add_definitions (-D_BSD_SOURCE)
  # Set some linker flags
  if(CMAKE_SYSTEM_NAME MATCHES Linux)
    set(CMAKE_SHARED_LINKER_FLAGS
      "-lc ${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS
      "-lc ${CMAKE_MODULE_LINKER_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS
      "-lc ${CMAKE_EXE_LINKER_FLAGS}")
  endif (CMAKE_SYSTEM_NAME MATCHES Linux)
  check_cxx_compiler_flag(-fvisibility=hidden HAVE_INTEL_VISIBILITY)
  set(HAVE_INTEL_VISIBILITY ${HAVE_INTEL_VISIBILITY} CACHE BOOL
    "Intel Compiler support for hidden visibility")
  if(HAVE_INTEL_VISIBILITY AND ENABLE_SYMBOL_VISIBILITY)
    add_definitions(-DHAVE_GCC_VISIBILITY)
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
    check_cxx_compiler_flag(-fvisibility-inlines-hidden HAVE_INTEL_VISIBILITY_INLINES)
    set(HAVE_INTEL_VISIBILITY_INLINES ${HAVE_INTEL_VISIBILITY_INLINES} CACHE BOOL
      "Intel Compiler support for -fvisibility-inlines-hidden")
    if(HAVE_INTEL_VISIBILITY_INLINES)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility-inlines-hidden")
    endif()
  endif()

  set(CMAKE_C_FLAGS_RELEASE "-O3 -ipo1 -DNDEBUG -Wl,-s")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ipo1 -DNDEBUG -Wl,-s")

  # define WITH_SSE2 to enable the SSE2 instruction set.
  # Available on Pentium 4, Athlon 64, and newer CPUs.
  if(CMAKE_SIZEOF_VOID_P EQUAL 8 OR WITH_SSE2)
    message(STATUS "Enabling the SSE2 instruction set")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
  endif(CMAKE_SIZEOF_VOID_P EQUAL 8 OR WITH_SSE2)
endif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel AND UNIX)

#Find if OpenMP support is enabled

option(ENABLE_OPENMP
    "Enable support for OpenMP compilation of forcefield code"
    OFF)
if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()

if(MSVC)
  message(STATUS "MSVC_VERSION =" ${MSVC_VERSION})
  if (${MSVC_VERSION} EQUAL 1600)
    # To be added: AND SP1  OR gcc>4.0
    add_definitions(-DHAVE_SHARED_POINTER)
    set(SHARED_POINTER true)
    message("SHARED_POINTER set") 
  endif()
endif()

# Some rpath handling for Linux and Mac
if(UNIX AND BUILD_SHARED)
  if(APPLE)
    set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR})
  else()
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif()
endif()

# this is adapted from KDE's FindKDE4Internal.cmake : default the build type to
# "release with debug info".
#
# We will define two other build types: Debug and Release.
# These names are case-insensitive i.e. you can do -DCMAKE_BUILD_TYPE=debug
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

#define various build types
if(CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2 -g")
  set (CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -g3 -fno-inline")
endif()

file(GLOB_RECURSE headers include/openbabel/*.h)
configure_file(${openbabel_SOURCE_DIR}/src/config.h.cmake
  ${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h)
install(FILES ${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h
  DESTINATION include/openbabel-2.0/openbabel)
if(NOT MSVC)
  configure_file(${openbabel_SOURCE_DIR}/openbabel-2.0.pc.cmake
    ${openbabel_BINARY_DIR}/openbabel-2.0.pc @ONLY)
  install(FILES ${openbabel_BINARY_DIR}/openbabel-2.0.pc
    DESTINATION ${LIB_INSTALL_DIR}/pkgconfig)
endif()

include_directories(${openbabel_BINARY_DIR}/include
  ${openbabel_SOURCE_DIR}/data
  ${openbabel_SOURCE_DIR}/include
)

#cjh
find_package(Eigen3)
if(EIGEN3_FOUND)
  add_definitions(-DHAVE_EIGEN -DHAVE_EIGEN3)
  include_directories(${EIGEN3_INCLUDE_DIR})
else()
  find_package(Eigen2) # find and setup Eigen2
  if(EIGEN2_FOUND)
    add_definitions (-DHAVE_EIGEN)
    include_directories(${EIGEN2_INCLUDE_DIR})
  endif()
endif()
#cjh

find_package(Cairo)
if(CAIRO_FOUND)
  message(STATUS "Cairo found. PNG output will be supported.")
else(CAIRO_FOUND)
  message(STATUS "Could NOT find Cairo. PNG output will NOT be supported.")
endif(CAIRO_FOUND)

add_subdirectory(include)
add_subdirectory(data)
add_subdirectory(doc)
add_subdirectory(src)
add_subdirectory(tools)

# Should we build the GUI? Default is yes on Linux or Windows, not Mac.
if (APPLE OR NOT BUILD_SHARED)
  option(BUILD_GUI "Build the GUI" OFF)
else ()
  option(BUILD_GUI "Build the GUI" ON)
endif()

if(BUILD_GUI)
  message(STATUS "Attempting to build the GUI")
  if(wxWidgets_FOUND)
    message(STATUS "   wxWidgets found => GUI will be built")
    add_subdirectory(src/GUI)
  else()
    message(STATUS "   wxWidgets not found => GUI will not be built")
  endif()
else()
  message(STATUS "GUI will not be built")
endif()

# Do we want to build the tests?
option(ENABLE_TESTS "Enable unit tests" ON)
if(ENABLE_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(test)
endif()

# Should the language bindings be regenereted?
option(RUN_SWIG "Generate language bindings with SWIG" OFF)

# Should all bindings be built?
option(ALL_BINDINGS "Build all languages bindings" OFF)

# Should Python bindings be built?
option(PYTHON_BINDINGS "Build Python bindings" OFF)

# Should Ruby bindings be built?
option(RUBY_BINDINGS "Build Ruby bindings" OFF)

# Should Perl bindings be built?
option(PERL_BINDINGS "Build Perl bindings" OFF)

# Should Java bindings be built?
option(JAVA_BINDINGS "Build Java bindings" OFF)

# Should CSharp bindings be built?
option(CSHARP_BINDINGS "Build Csharp bindings" OFF)

if(ALL_BINDINGS)
  set(PYTHON_BINDINGS ON)
  set(RUBY_BINDINGS ON)
  set(PERL_BINDINGS ON)
  set(JAVA_BINDINGS ON)
  set(CSHARP_BINDINGS ON)
endif()

add_subdirectory(scripts)

# A minimal build with just two formats (SMI, SDF) useful for developers
option(MINIMAL_BUILD "Disable full build" OFF)

# Create the Config and ConfigVersion files in the build directory, useful to
# build against an uninstalled OpenBabel build. Simply set OpenBabel2_DIR to
# point at the build directory and call find_package(OpenBabel2)
set(OpenBabel2_INCLUDE_DIRS
  "${openbabel_SOURCE_DIR}/include"
  "${openbabel_BINARY_DIR}/include")
configure_file("${openbabel_SOURCE_DIR}/OpenBabel2Config.cmake.in"
  "${openbabel_BINARY_DIR}/OpenBabel2Config.cmake" @ONLY)
configure_file("${openbabel_SOURCE_DIR}/OpenBabel2ConfigVersion.cmake.in"
  "${openbabel_BINARY_DIR}/OpenBabel2ConfigVersion.cmake" @ONLY)

# Now to configure the installed config file.
set(OB_CONFIG_DIR "${LIB_INSTALL_DIR}/cmake/openbabel2")
set(OpenBabel2_INCLUDE_DIRS "\${OpenBabel2_INSTALL_PREFIX}/${OB_INCLUDE_DIRS}")
set(OB_EXPORTS_FILE
  "\${OpenBabel2_INSTALL_PREFIX}/${OB_CONFIG_DIR}/OpenBabel2_EXPORTS.cmake")
set(REL_REF)
if(NOT WIN32)
  set(REL_REF "/../../..")
endif()
configure_file("${openbabel_SOURCE_DIR}/OpenBabel2Config.cmake.in"
  "${openbabel_BINARY_DIR}/OpenBabel2ConfigInstall.cmake" @ONLY)

install(EXPORT openbabel-targets
  DESTINATION "${OB_CONFIG_DIR}"
  FILE OpenBabel2_EXPORTS.cmake
  COMPONENT Development)
install(FILES "${openbabel_BINARY_DIR}/OpenBabel2ConfigInstall.cmake"
  DESTINATION ${OB_CONFIG_DIR}
  COMPONENT Development
  RENAME OpenBabel2Config.cmake)
install(FILES "${openbabel_BINARY_DIR}/OpenBabel2ConfigVersion.cmake"
  DESTINATION ${OB_CONFIG_DIR}
  COMPONENT Development)
