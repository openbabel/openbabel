# Please ensure that any changes remain compliant with 3.1.
if(NOT EMBED_OPENBABEL)
  cmake_minimum_required(VERSION 3.1)
endif()

project(openbabel)
set(CMAKE_MODULE_PATH ${openbabel_SOURCE_DIR}/cmake/modules)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  if(POLICY CMP0042)
    cmake_policy(SET CMP0042 OLD)
  endif()
endif()

include (CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Multi-threaded applications may experience race conditions.")
endif()

#include (MacroEnsureVersion)

# Version numbering - should be bumped for each release
# Note that for "beta" releases, we should start at x.90.0 -- we've
# had too many 1.100.1 releases. :-)
set(BABEL_MAJ_VER   3)
set(BABEL_MIN_VER   0)
set(BABEL_PATCH_VER 0)

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
if(NOT DEFINED OB_INCLUDE_DIRS)
  set(OB_INCLUDE_DIRS "include/openbabel${BABEL_MAJ_VER}")
endif()
set(OB_EXPORTS_FILE "${openbabel_BINARY_DIR}/OpenBabel3_EXPORTS.cmake")
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

#GLIBC_2.4
option(GLIBC_24_COMPATIBLE "Build project compatible with GLIBC 2.4" OFF)
if(GLIBC_24_COMPATIBLE AND UNIX)
  set(WRAP_MEMCPY_SRC memcpy.c)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--wrap=memcpy")
endif()

option(BUILD_SHARED "enable shared build support" ON)
option(BUILD_MIXED "enable linux relocatable binaries support" OFF)
option(WITH_STATIC_INCHI "enable static inchi and related plugins" OFF)
option(WITH_STATIC_LIBXML "enable static libxml and xml file format plugins" OFF)
# Static building
if(BUILD_SHARED)
  set(BUILD_TYPE SHARED)
  set(PLUGIN_TYPE MODULE)
  add_definitions(-DUSING_DYNAMIC_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
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
    set(CMAKE_C_FLAGS_DEBUG "/MTd")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "/INCREMENTAL:NO /NODEFAULTLIB:MSVCRT")
    set(CMAKE_CXX_FLAGS_RELEASE	"/MT /O2 /Ob2 /D NDEBUG")
    set(CMAKE_C_FLAGS_RELEASE	"/MT /O2 /Ob2 /D NDEBUG")
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
  set(INCHI_LIBRARY inchi)
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
  include_directories(${openbabel_SOURCE_DIR}/windows-msvc/include)
endif()

if(MSVC AND OB_USE_PREBUILT_BINARIES)
  include_directories(${XDR_INCLUDE_DIR})
  # The following line is required for check_include_file(rpc/xdr.h ...) below
  set(CMAKE_REQUIRED_INCLUDES ${XDR_INCLUDE_DIR})
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

# Cygwin may be missing an XDR function: https://www.gnu.org/software/gnulib/manual/html_node/xdrstdio_005fcreate.html
if(CYGWIN)
  set(HAVE_RPC_XDR_H FALSE)
  message(STATUS "Disabling rpc/xdr.h for Cygwin")
endif()

# Get the GCC version - from KDE4 cmake files
if(CMAKE_COMPILER_IS_GNUCXX)
  if(NOT(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.0.0))
    set(GCC_IS_NEWER_THAN_4_0 TRUE)
  else()
    set(GCC_IS_NEWER_THAN_4_0 FALSE)
  endif()
  if(NOT(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.1.0))
    set(GCC_IS_NEWER_THAN_4_1 TRUE)
  else()
    set(GCC_IS_NEWER_THAN_4_1 FALSE)
  endif()
  if(NOT(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.2.0))
    set(GCC_IS_NEWER_THAN_4_2 TRUE)
  else()
    set(GCC_IS_NEWER_THAN_4_2 FALSE)
  endif()
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

message(STATUS "Checking available shared pointer...")
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <memory>
                           int main(){ std::shared_ptr<int> p; return 0; }"
                           HAVE_STD_SHARED_PTR)
check_cxx_source_compiles("#include <memory>
                           int main(){ std::tr1::shared_ptr<int> p; return 0; }"
                           HAVE_STD_TR1_SHARED_PTR)
check_cxx_source_compiles("#include <tr1/memory>
                           int main(){ std::tr1::shared_ptr<int> p; return 0; }"
                           HAVE_STD_TR1_SHARED_PTR_IN_TR1_HEADER)
if(HAVE_STD_SHARED_PTR)
  set(OB_SHARED_PTR_IMPLEMENTATION "std::shared_ptr")
  set(OB_SHARED_PTR_HEADER "memory")
  set(SHARED_POINTER true)
  add_definitions(-DHAVE_SHARED_POINTER)
elseif(HAVE_STD_TR1_SHARED_PTR)
  set(OB_SHARED_PTR_IMPLEMENTATION "std::tr1::shared_ptr")
  set(OB_SHARED_PTR_HEADER "memory")
  set(SHARED_POINTER true)
  add_definitions(-DHAVE_SHARED_POINTER)
elseif(HAVE_STD_TR1_SHARED_PTR_IN_TR1_HEADER)
  set(OB_SHARED_PTR_IMPLEMENTATION "std::tr1::shared_ptr")
  set(OB_SHARED_PTR_HEADER "tr1/memory")
  set(SHARED_POINTER true)
  add_definitions(-DHAVE_SHARED_POINTER)
else()
  find_package(Boost)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
    add_definitions(-DUSE_BOOST -DHAVE_SHARED_POINTER)
    set(SHARED_POINTER true)
  else()
    message(FATAL_ERROR "Boost is required to compile OpenBabel with GCC 3.x")
  endif()
  set(OB_SHARED_PTR_IMPLEMENTATION "boost::shared_ptr")
  set(OB_SHARED_PTR_HEADER "boost/memory")
endif()
message(STATUS "Use ${OB_SHARED_PTR_IMPLEMENTATION} in ${OB_SHARED_PTR_HEADER}")

find_package(Boost COMPONENTS program_options iostreams)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
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
  # TODO: Remove this? Should be handled by -march=native and only enabled if OPTIMIZE_NATIVE is ON
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

# Three build types are available: Release, Debug, RelWithDebInfo.
# We default to Release.
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release)
endif()

#define various build types
if(CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2 -g")
  set (CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -g3 -fno-inline")
endif()

file(GLOB_RECURSE headers include/openbabel/*.h)
configure_file(${openbabel_SOURCE_DIR}/src/config.h.cmake
  ${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h)
install(FILES ${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h
  DESTINATION ${OB_INCLUDE_DIRS}/openbabel)
if(NOT MSVC)
  configure_file(${openbabel_SOURCE_DIR}/openbabel-3.pc.cmake
    ${openbabel_BINARY_DIR}/openbabel-3.pc @ONLY)
  install(FILES ${openbabel_BINARY_DIR}/openbabel-3.pc
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


find_package(Boost)
if(Boost_FOUND AND BUILD_SHARED)
    include_directories(${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
    option(WITH_MAEPARSER "Build Maestro support" ON)
else()
    option(WITH_MAEPARSER "Build Maestro support" OFF)
endif()


if(WITH_MAEPARSER)
    set (CMAKE_CXX_STANDARD 11)

    set(MAEPARSER_VERSION 1.1)
    if(NOT EXISTS "${openbabel_SOURCE_DIR}/external/maeparser-${MAEPARSER_VERSION}")
        file(DOWNLOAD "https://github.com/schrodinger/maeparser/archive/v${MAEPARSER_VERSION}.tar.gz"
            "${openbabel_SOURCE_DIR}/external/maeparser-${MAEPARSER_VERSION}.tar.gz" STATUS status)
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
            ${openbabel_SOURCE_DIR}/external/maeparser-${MAEPARSER_VERSION}.tar.gz
            WORKING_DIRECTORY ${openbabel_SOURCE_DIR}/external)
        message(STATUS "Downloaded MaeParser to ${openbabel_SOURCE_DIR}/external/maeparser-${MAEPARSER_VERSION}.")
    endif()
    set(MAEPARSER_DIR ${openbabel_SOURCE_DIR}/external/maeparser-${MAEPARSER_VERSION})
    add_subdirectory(${MAEPARSER_DIR})
    include_directories(${MAEPARSER_DIR})
    set(libs ${libs} maeparser)

    message(STATUS "Maestro formats will be supported. Using MaeParser ${MAEPARSER_VERSION} at ${MAEPARSER_INCLUDE_DIRS}")
else()
    message(STATUS "Maestro formats will NOT be supported. Please install Boost to enable Maestro formats.")
endif()

option(WITH_JSON "Build JSON support" ON)
option(OPENBABEL_USE_SYSTEM_RAPIDJSON "Use the system rapidjson if available." ON)
if(WITH_JSON)
  set(RAPIDJSON_VERSION_MIN 1.1.0)
  # First try find system rapidjson
  if(OPENBABEL_USE_SYSTEM_RAPIDJSON)
    find_package(RapidJSON ${RAPIDJSON_VERSION_MIN})
  endif()
  # Otherwise download rapidjson to openbabel source directory
  if (NOT RAPIDJSON_FOUND OR RAPIDJSON_VERSION VERSION_LESS ${RAPIDJSON_VERSION_MIN})
    set(RAPIDJSON_VERSION 1.1.0)
    if(NOT EXISTS "${openbabel_SOURCE_DIR}/external/rapidjson-${RAPIDJSON_VERSION}")
      file(DOWNLOAD "https://github.com/Tencent/rapidjson/archive/v${RAPIDJSON_VERSION}.tar.gz"
        "${openbabel_SOURCE_DIR}/external/rapidjson-${RAPIDJSON_VERSION}.tar.gz" STATUS status)
      execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
        ${openbabel_SOURCE_DIR}/external/rapidjson-${RAPIDJSON_VERSION}.tar.gz
        WORKING_DIRECTORY ${openbabel_SOURCE_DIR}/external)
      message(STATUS "Downloaded RapidJSON to ${openbabel_SOURCE_DIR}/external/rapidjson-${RAPIDJSON_VERSION}.")
    endif()
    set(RAPIDJSON_FOUND TRUE)
    set(RAPIDJSON_INCLUDE_DIRS ${openbabel_SOURCE_DIR}/external/rapidjson-${RAPIDJSON_VERSION}/include)
  endif()
  if(RAPIDJSON_FOUND)
    message(STATUS "JSON formats will be supported. Using RapidJSON ${RAPIDJSON_VERSION} at ${RAPIDJSON_INCLUDE_DIRS}")
    add_definitions(-DHAVE_RADPIJSON)
    include_directories(${RAPIDJSON_INCLUDE_DIRS})
  endif()
endif()

find_package(Cairo)
if(CAIRO_FOUND)
  message(STATUS "Cairo found. PNG output will be supported.")
else(CAIRO_FOUND)
  message(STATUS "Could NOT find Cairo. PNG output will NOT be supported.")
endif(CAIRO_FOUND)

# Should we enable optimizations for the native CPU architecture?
# (this will speed up JSON handling, similarity calculations and maybe those involving Eigen)
option(OPTIMIZE_NATIVE "Optimize for native CPU architecture. Turn off if compiling for distribution/reuse on other machines." OFF)
if(OPTIMIZE_NATIVE)
  include(OptimizeNative)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NATIVE_FLAGS}")
  # Tell rapidjson which extensions are available
  if(HAVE_SSE2)
    add_definitions(-DRAPIDJSON_SSE2)
  endif()
  if(HAVE_SSE4_2)
    add_definitions(-DRAPIDJSON_SSE42)
  endif()
  if(HAVE_NEON)
    add_definitions(-DRAPIDJSON_NEON)
  endif()
  add_definitions(-DOPTIMIZE_NATIVE)
  message(STATUS "Optimizing code for this machine's CPU architecture. Use -DOPTIMIZE_NATIVE=OFF if compiling for distribution/reuse on other machines.")
endif()

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

# Build bindings only
option(BINDINGS_ONLY "Build bindings only" OFF)

# Point to library if building bindings only
set(BABEL_SYSTEM_LIBRARY ${BABEL_LIBRARY}
    CACHE PATH "Point to openbabel library if building bindings only")

# Should all bindings be built?
option(ALL_BINDINGS "Build all languages bindings" OFF)

# Should PHP bindings be built
option(PHP_BINDINGS "Build PHP bindings" OFF)

# Should Python bindings be built?
option(PYTHON_BINDINGS "Build Python bindings" OFF)

# Should Ruby bindings be built?
option(RUBY_BINDINGS "Build Ruby bindings" OFF)

# Should Perl bindings be built?
option(PERL_BINDINGS "Build Perl bindings" OFF)

# Should Java bindings be built?
option(JAVA_BINDINGS "Build Java bindings" OFF)

# Should R bindings be built?
option(R_BINDINGS "Build R bindings" OFF)

# Should CSharp bindings be built?
option(CSHARP_BINDINGS "Build Csharp bindings" OFF)

if(ALL_BINDINGS)
  set(PYTHON_BINDINGS ON)
  set(RUBY_BINDINGS ON)
  set(PHP_BINDINGS ON)
  set(PERL_BINDINGS ON)
  set(JAVA_BINDINGS ON)
  set(CSHARP_BINDINGS ON)
  set(R_BINDINGS ON)
endif()

add_subdirectory(scripts)

# A minimal build with just two formats (SMI, SDF) useful for developers
option(MINIMAL_BUILD "Disable full build" OFF)

# Create the Config and ConfigVersion files in the build directory, useful to
# build against an uninstalled OpenBabel build. Simply set OpenBabel3_DIR to
# point at the build directory and call find_package(OpenBabel3)
set(OpenBabel3_INCLUDE_DIRS
  "${openbabel_SOURCE_DIR}/include"
  "${openbabel_BINARY_DIR}/include")
configure_file("${openbabel_SOURCE_DIR}/OpenBabel3Config.cmake.in"
  "${openbabel_BINARY_DIR}/OpenBabel3Config.cmake" @ONLY)
configure_file("${openbabel_SOURCE_DIR}/OpenBabel3ConfigVersion.cmake.in"
  "${openbabel_BINARY_DIR}/OpenBabel3ConfigVersion.cmake" @ONLY)

# Now to configure the installed config file.
set(OB_CONFIG_DIR "${LIB_INSTALL_DIR}/cmake/openbabel3")
set(OpenBabel3_INCLUDE_DIRS "\${OpenBabel3_INSTALL_PREFIX}/${OB_INCLUDE_DIRS}")
set(OB_EXPORTS_FILE
  "\${OpenBabel3_INSTALL_PREFIX}/${OB_CONFIG_DIR}/OpenBabel3_EXPORTS.cmake")
set(REL_REF)
if(NOT WIN32)
  set(REL_REF "/../../..")
endif()
configure_file("${openbabel_SOURCE_DIR}/OpenBabel3Config.cmake.in"
  "${openbabel_BINARY_DIR}/OpenBabel3ConfigInstall.cmake" @ONLY)

install(EXPORT openbabel-targets
  DESTINATION "${OB_CONFIG_DIR}"
  FILE OpenBabel3_EXPORTS.cmake
  COMPONENT Development)
install(FILES "${openbabel_BINARY_DIR}/OpenBabel3ConfigInstall.cmake"
  DESTINATION ${OB_CONFIG_DIR}
  COMPONENT Development
  RENAME OpenBabel3Config.cmake)
install(FILES "${openbabel_BINARY_DIR}/OpenBabel3ConfigVersion.cmake"
  DESTINATION ${OB_CONFIG_DIR}
  COMPONENT Development)

# from https://cmake.org/Wiki/RecipeAddUninstallTarget
########### Add uninstall target ###############
if (NOT DEFINED EMBED_OPENBABEL)
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake")
endif()
