# Usage:
# SET_SOURCE_FILES_PROPERTIES(test.py PROPERTIES PYTHONPATH
#   "${LIBRARY_OUTPUT_PATH}:${VTK_DIR}")
# ADD_PYTHON_TEST(PYTHON-TEST test.py)
#
#  Copyright (c) 2006-2010 Mathieu Malaterre <mathieu.malaterre@gmail.com>
#
#  Redistribution and use is allowed according to the terms of the New
#  BSD license.
#  For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

# Set CMake Policy regarding the use of the @VARNAME@ below
if(POLICY CMP0053)
  cmake_policy(PUSH) # Note the corresponding POP at the end of the file
  cmake_policy(SET CMP0053 NEW)
endif()

# Need python interpreter:
FIND_PACKAGE(PythonInterp)
MARK_AS_ADVANCED(PYTHON_EXECUTABLE)

# Make sure we handle systems w/o python (e.g. chroot)
if(NOT PYTHON_EXECUTABLE)
  if(POLICY CMP0053)
    cmake_policy(POP)
  endif()
  return()
endif(NOT PYTHON_EXECUTABLE)

MACRO(ADD_PYTHON_TEST TESTNAME FILENAME)
  GET_SOURCE_FILE_PROPERTY(loc ${FILENAME} LOCATION)
  GET_SOURCE_FILE_PROPERTY(pyenv ${FILENAME} PYTHONPATH)
  GET_SOURCE_FILE_PROPERTY(ob_libdir ${FILENAME} BABEL_LIBDIR)
  GET_SOURCE_FILE_PROPERTY(ob_datadir ${FILENAME} BABEL_DATADIR)
  STRING(REGEX REPLACE ";" " " wo_semicolumn "${ARGN}")
  ADD_TEST(NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${loc} ${wo_semicolumn})
  SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES ENVIRONMENT "PYTHONPATH=${pyenv};LD_LIBRARY_PATH=${pyenv}:$ENV{LD_LIBRARY_PATH};BABEL_LIBDIR=${ob_libdir};BABEL_DATADIR=${ob_datadir}")
ENDMACRO(ADD_PYTHON_TEST)

# Byte compile recursively a directory (DIRNAME)
MACRO(ADD_PYTHON_COMPILEALL_TEST DIRNAME)
  # First get the path:
  GET_FILENAME_COMPONENT(temp_path "${PYTHON_LIBRARIES}" PATH)
  # Find the python script:
  GET_FILENAME_COMPONENT(PYTHON_COMPILE_ALL_PY "${temp_path}/../compileall.py" ABSOLUTE)
  # add test, use DIRNAME to create uniq name for the test:
  ADD_TEST(COMPILE_ALL-${DIRNAME} ${PYTHON_EXECUTABLE} "${PYTHON_COMPILE_ALL_PY}" -q ${DIRNAME})
ENDMACRO(ADD_PYTHON_COMPILEALL_TEST)

if(POLICY CMP0053)
  cmake_policy(POP)
endif()
