cmake_minimum_required(VERSION 2.6.0)

include(FindSWIG)
# This module finds an installed SWIG. It sets the following variables:
# SWIG_FOUND - set to true if SWIG is found
# SWIG_DIR - the directory where swig is installed
# SWIG_EXECUTABLE - the path to the swig executable
# SWIG_VERSION   - the version number of the swig executable

find_package(SWIG)
if (NOT SWIG_FOUND)
  message(STATUS "SWIG NOT found")
  return()
endif (NOT SWIG_FOUND)
message(STATUS "SWIG version is ${SWIG_VERSION}")
if (NOT ${SWIG_VERSION} EQUAL "1.3.40")
  message(STATUS "SWIG version 1.3.40 exactly is required")
  return()
endif (NOT ${SWIG_VERSION} EQUAL "1.3.40")


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

add_custom_command(OUTPUT openbabel-python.cpp openbabel.py
    COMMAND ${SWIG_EXECUTABLE} -python -c++ -small -O -templatereduce -naturalvar -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o openbabel-python.cpp ${openbabel_SOURCE_DIR}/scripts/openbabel-python.i
    MAIN_DEPENDENCY openbabel-python.i
)

if(NOT WIN32)
    add_custom_target(_openbabel ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${openbabel_SOURCE_DIR}/scripts/python ${openbabel_BINARY_DIR}/scripts/python
        COMMAND ${PYTHON_EXECUTABLE} ${openbabel_BINARY_DIR}/scripts/python/setup.py build
        DEPENDS openbabel-python.cpp
        WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts
    )
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
        COMMAND ${CMAKE_COMMAND} -E remove ob.py
        VERBATIM)

    if (${PYTHON_PREFIX})
        install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/python/setup.py install --prefix=${PYTHON_PREFIX} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/scripts)")
    else (${PYTHON_PREFIX})
        if (${CMAKE_INSTALL_PREFIX} EQUAL "/usr/local")
            install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/python/setup.py install WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/scripts)")
        else (${CMAKE_INSTALL_PREFIX} EQUAL "/usr/local")
            install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/python/setup.py install --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/scripts)")
        endif (${CMAKE_INSTALL_PREFIX} EQUAL "/usr/local")
    endif (${PYTHON_PREFIX})

else(NOT WIN32)
    add_library(_openbabel MODULE openbabel-python.cpp)
    target_link_libraries(_openbabel ${PYTHON_LIBRARIES} ${BABEL_LIBRARY})
    set_target_properties(_openbabel PROPERTIES
	LIBRARY_OUTPUT_DIRECTORY ${openbabel_BINARY_DIR}
        OUTPUT_NAME _openbabel
        PREFIX ""
        SUFFIX .pyd )
    add_custom_command(TARGET _openbabel POST_BUILD 
      COMMAND ${CMAKE_COMMAND} -E copy openbabel.py ${openbabel_BINARY_DIR}/${CMAKE_CFG_INTDIR}
      )
endif(NOT WIN32)

