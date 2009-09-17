cmake_minimum_required(VERSION 2.6.0)
# This module finds an installed SWIG. It sets the following variables:
# SWIG_FOUND - set to true if SWIG is found
# SWIG_DIR - the directory where swig is installed
# SWIG_EXECUTABLE - the path to the swig executable
# SWIG_VERSION   - the version number of the swig executable
include(FindSWIG)
# Defines the following macros:
# SWIG_ADD_MODULE(name language [ files ])
#     - Define swig module with given name and specified language
# SWIG_LINK_LIBRARIES(name [ libraries ])
#     - Link libraries to swig module
include(UseSWIG)

find_package(SWIG)
if (NOT SWIG_FOUND)
  message(STATUS "SWIG NOT found")
  return()
endif (NOT SWIG_FOUND)

###################
# Python
###################

find_package(PythonLibs)
if (NOT PYTHONLIBS_FOUND)
  message(STATUS "Python libraries NOT found")
  return()
endif (NOT PYTHONLIBS_FOUND)

find_package(PythonInterp)
if (NOT PYTHONINTERP_FOUND)
  message(STATUS "Python interpreter NOT found")
  return()
endif (NOT PYTHONINTERP_FOUND)

include_directories(
    ${PYTHON_INCLUDE_PATH}
)

# python bindings
set_source_files_properties(openbabel-python.i PROPERTIES CPLUSPLUS ON)
set(CMAKE_SWIG_FLAGS -small -O -python -templatereduce -naturalvar)
swig_add_module(_openbabel c++ openbabel-python.i)
swig_link_libraries(_openbabel openbabel)

if(NOT WIN32)
add_custom_command(TARGET _openbabel POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E echo "import sys" > ob.py
    COMMAND ${CMAKE_COMMAND} -E echo "if sys.platform.find('linux'\) != -1:" >> ob.py
    COMMAND ${CMAKE_COMMAND} -E echo "  try:" >> ob.py
    COMMAND ${CMAKE_COMMAND} -E echo "    import dl" >> ob.py
    COMMAND ${CMAKE_COMMAND} -E echo "  except ImportError:" >> ob.py
    COMMAND ${CMAKE_COMMAND} -E echo "    import DLFCN as dl" >> ob.py
    COMMAND ${CMAKE_COMMAND} -E echo "  sys.setdlopenflags(sys.getdlopenflags() | dl.RTLD_GLOBAL)" >> ob.py
    COMMAND cat openbabel.py >> ob.py
    COMMAND ${CMAKE_COMMAND} -E copy ob.py openbabel.py
    VERBATIM)
endif(NOT WIN32)

# find site-packages dir...
if(NOT PYTHON_PREFIX)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c 
    "from sys import stdout; from distutils import sysconfig; stdout.write(sysconfig.get_python_lib())"
    OUTPUT_VARIABLE PYTHON_LIB_PATH)
else(NOT PYTHON_PREFIX)
  set(PYTHON_LIB_PATH ${PYTHON_PREFIX})
endif(NOT PYTHON_PREFIX)

message(STATUS "Python site-packages directory: ${PYTHON_LIB_PATH}")

install(TARGETS _openbabel DESTINATION ${PYTHON_LIB_PATH})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/openbabel.py DESTINATION ${PYTHON_LIB_PATH})


