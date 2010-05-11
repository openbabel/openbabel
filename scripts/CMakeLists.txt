cmake_minimum_required(VERSION 2.6.0)

if (RUN_SWIG)
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
  if (NOT ${SWIG_VERSION} STREQUAL "1.3.40")
    message(STATUS "SWIG version 1.3.40 exactly is required")
    return()
  endif (NOT ${SWIG_VERSION} STREQUAL "1.3.40")
endif (RUN_SWIG)

###################
# Python
###################

if (PYTHON_BINDINGS)
  if (EXISTS ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp OR
      RUN_SWIG)
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

    if (RUN_SWIG)
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp ${openbabel_SOURCE_DIR}/scripts/python/openbabel.py
          COMMAND ${SWIG_EXECUTABLE} -python -c++ -small -O -templatereduce -naturalvar -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp -outdir ${openbabel_SOURCE_DIR}/scripts/python ${openbabel_SOURCE_DIR}/scripts/openbabel-python.i
          MAIN_DEPENDENCY openbabel-python.i
      )
    endif (RUN_SWIG)

    if(NOT WIN32)
        add_custom_target(_openbabel ALL
            COMMAND ${PYTHON_EXECUTABLE} ${openbabel_SOURCE_DIR}/scripts/python/setup.py build --build-lib=${openbabel_BINARY_DIR}/scripts/pybuild
            DEPENDS ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp
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
            COMMAND cat ${openbabel_SOURCE_DIR}/scripts/python/openbabel.py >> ob.py
            COMMAND ${CMAKE_COMMAND} -E copy ob.py ${openbabel_SOURCE_DIR}/scripts/python/openbabel.py
            COMMAND ${CMAKE_COMMAND} -E remove ob.py
            VERBATIM)

        if (${PYTHON_PREFIX})
            install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${openbabel_SOURCE_DIR}/scripts/python/setup.py install --prefix=${PYTHON_PREFIX} WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts)")
        else (${PYTHON_PREFIX})
            if (${CMAKE_INSTALL_PREFIX} EQUAL "/usr/local")
                install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${openbabel_SOURCE_DIR}/scripts/python/setup.py install WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts)")
            else (${CMAKE_INSTALL_PREFIX} EQUAL "/usr/local")
                install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${openbabel_SOURCE_DIR}/scripts/python/setup.py install --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts)")
            endif (${CMAKE_INSTALL_PREFIX} EQUAL "/usr/local")
        endif (${PYTHON_PREFIX})

    else(NOT WIN32)
        add_library(bindings_python MODULE ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp)
        target_link_libraries(bindings_python ${PYTHON_LIBRARIES} ${BABEL_LIBRARY})
        set_target_properties(bindings_python PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${openbabel_BINARY_DIR}
            OUTPUT_NAME _openbabel
            PREFIX ""
            SUFFIX .pyd )
        add_custom_command(TARGET bindings_python POST_BUILD 
          COMMAND ${CMAKE_COMMAND} -E copy ${openbabel_SOURCE_DIR}/scripts/python/openbabel.py ${openbabel_BINARY_DIR}/${CMAKE_CFG_INTDIR}
          )
    endif(NOT WIN32)
  else (EXISTS ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp OR
      RUN_SWIG)
    message(STATUS "Warning: Python bindings NOT found. Generate using -DRUN_SWIG.")
  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp OR
      RUN_SWIG)
endif (PYTHON_BINDINGS)

###################
# Java
###################

if (JAVA_BINDINGS)
  if (EXISTS ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp OR
      RUN_SWIG)

    find_package(JNI)
    if (NOT JNI_INCLUDE_DIRS)
      message(STATUS "Java JNI include files NOT found. Java bindings will NOT be compiled.")
      return()
    endif (NOT JNI_INCLUDE_DIRS)
    include_directories(
        ${JNI_INCLUDE_DIRS}
    )

    if (RUN_SWIG)
      find_package(Java)
      if (NOT JAVA_COMPILE)
        message(STATUS "Java compiler NOT found. Java bindings will NOT be generated.")
        return()
      endif (NOT JAVA_COMPILE)
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp # ${openbabel_SOURCE_DIR}/scripts/python/openbabel.py
          COMMAND ${CMAKE_COMMAND} -E make_directory ${openbabel_SOURCE_DIR}/scripts/java/org/openbabel
          COMMAND ${SWIG_EXECUTABLE} -java -package org.openbabel -c++ -small -O -templatereduce -naturalvar -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp -outdir ${openbabel_SOURCE_DIR}/scripts/java/org/openbabel ${openbabel_SOURCE_DIR}/scripts/openbabel-java.i
          COMMAND ${JAVA_COMPILE} ${openbabel_SOURCE_DIR}/scripts/java/org/openbabel/*.java
          COMMAND ${CMAKE_COMMAND} -E chdir ${openbabel_SOURCE_DIR}/scripts/java ${JAVA_ARCHIVE} cf openbabel.jar org
          COMMAND ${CMAKE_COMMAND} -E remove_directory ${openbabel_SOURCE_DIR}/scripts/java/org
          MAIN_DEPENDENCY openbabel-java.i
      )
    endif (RUN_SWIG)

    add_library(bindings_java MODULE ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp)
    target_link_libraries(bindings_java ${JNI_LIBRARIES} ${BABEL_LIBRARY})
    set_target_properties(bindings_java PROPERTIES
        OUTPUT_NAME openbabel_java )
    if(MSVC)
        set_target_properties(bindings_java PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${openbabel_BINARY_DIR})
    endif(MSVC)

  else (EXISTS ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp OR
      RUN_SWIG)
    message(STATUS "Warning: Java bindings NOT found. Generate using -DRUN_SWIG.")
  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp OR
      RUN_SWIG)
endif (JAVA_BINDINGS)
