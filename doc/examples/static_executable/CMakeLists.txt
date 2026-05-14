#
# This script can be used to create static executables linking to the static
# OpenBabel3 library.
#
# This script requires OpenBabel to be build and installed. For example:
#
#  cd openbabel
#  mkdir build
#  cd build
#  cmake -DBUILD_SHARED=OFF -DCMAKE_INSTALL_PREFIX=/home/me/some/path ..
#  make
#  make install
#
# To compile your static executable:
#
#  cd myproject
#  mkdir build
#  cd build
#  cmake -DOpenBabel3_DIR=/home/me/some/path/lib/openbabel ..
#  make
#
# All plugins are inside the static libopenbabel.a but the symbols for the
# plugin classes have to be undefined. Plugins can be disabled by removing
# the class names from the format_classes, descriptor_classes, ... lists below.
#

# This line is required for cmake backwards compatibility.
cmake_minimum_required(VERSION 3.10)

# Name of your project
project(myproject)
# Create a list of source files (easier to maintain)
set(sources myexe.cpp)
# Set the name for the executable
set(executable_target myexe)

################################################################################

#
# Set compile flags for various compilers.
#
if(MSVC)
  # Set cl flags for static compiling
  set(CMAKE_CXX_FLAGS_DEBUG "/MTd")
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "/INCREMENTAL:NO /NODEFAULTLIB:MSVCRT")
  set(CMAKE_CXX_FLAGS_RELEASE "/MT /O2 /Ob2 /D NDEBUG")
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/INCREMENTAL:NO /NODEFAULTLIB:MSVCRT")
  # Note: static libraries are specified when running cmake
else()
  # Use -static flag to create static executable
  set(CMAKE_CXX_FLAGS "-static ${CMAKE_CXX_FLAGS}")
  # Make sure we find static libraries
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

# Set the path containing OpenBabel3Config.cmake, needed for find_package below.
find_path(OpenBabel3_DIR OpenBabel3Config.cmake PATHS
          ${OpenBabel3_DIR}
          "/usr/lib/openbabel"
         "/usr/local/lib/openbabel")

#
# Find and setup OpenBabel3.
#
find_package(OpenBabel3 REQUIRED)
include_directories(${OpenBabel3_INCLUDE_DIRS})

# Dependencies
find_package(LibXml2)

# The executable
add_executable(${executable_target} ${sources})
# Link against imported openbabel target
target_link_libraries(${executable_target} openbabel ${LIBXML2_LIBRARIES})
# Prevent -Wl,-Bdynamic from being added to the end of the link line.
set_target_properties(${executable_target} PROPERTIES
                      LINK_SEARCH_END_STATIC TRUE)
install(TARGETS ${executable_target} DESTINATION bin)
