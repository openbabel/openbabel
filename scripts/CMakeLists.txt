cmake_minimum_required(VERSION 2.6.0)
# Library versioning (used in Mac Python bindings)x
set(SOVERSION 4)

if (RUN_SWIG)
  include(FindSWIG)
  # This module finds an installed SWIG. It sets the following variables:
  # SWIG_FOUND - set to true if SWIG is found
  # SWIG_DIR - the directory where swig is installed
  # SWIG_EXECUTABLE - the path to the swig executable
  # SWIG_VERSION   - the version number of the swig executable

  find_package(SWIG 2.0)
  if (NOT SWIG_FOUND)
    message(STATUS "SWIG 2.0 NOT found")
    return()
  endif (NOT SWIG_FOUND)
endif (RUN_SWIG)

# The following needs to be kept in synch with ..\CMakeLists.txt
if(EIGEN3_FOUND)
  set (eigen_define "-DHAVE_EIGEN" "-DHAVE_EIGEN3")
  include_directories(${EIGEN3_INCLUDE_DIR})
  set (CPPFLAGS_VALUE "-I${EIGEN3_INCLUDE_DIR}")
  set (CPPFLAGS_TEXT "CPPFLAGS=${CPPFLAGS_VALUE}")
else()
  if(EIGEN2_FOUND)
    set (eigen_define "-DHAVE_EIGEN")
    include_directories(${EIGEN2_INCLUDE_DIR})
    set (CPPFLAGS_VALUE "-I${EIGEN2_INCLUDE_DIR}")
    set (CPPFLAGS_TEXT "CPPFLAGS=${CPPFLAGS_VALUE}")
  endif()
endif()

###################
# Python
###################

if (PYTHON_BINDINGS)
  if (EXISTS ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp OR RUN_SWIG)
    find_package(PythonInterp)
    if (NOT PYTHONINTERP_FOUND)
      message(STATUS "Python interpreter NOT found")
    endif (NOT PYTHONINTERP_FOUND)

    if (PYTHONINTERP_FOUND AND NOT Python_ADDITIONAL_VERSIONS)
      set(Python_ADDITIONAL_VERSIONS "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
    endif (PYTHONINTERP_FOUND AND NOT Python_ADDITIONAL_VERSIONS)

    find_package(PythonLibs)
    if (NOT PYTHONLIBS_FOUND)
      message(STATUS "Python libraries NOT found")
    endif (NOT PYTHONLIBS_FOUND)

    if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
      set(DO_PYTHON_BINDINGS ON BOOL)
      message(STATUS "Python bindings will be compiled")
    endif(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)

  else (EXISTS ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp OR RUN_SWIG)
    message(STATUS "Warning: Python bindings NOT found. Generate using -DRUN_SWIG=ON.")

  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp OR RUN_SWIG)
endif (PYTHON_BINDINGS)

if (DO_PYTHON_BINDINGS)

    include_directories(
        ${PYTHON_INCLUDE_PATH}
    )

    if(RUN_SWIG)
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/python/openbabel/openbabel-python.cpp ${openbabel_SOURCE_DIR}/scripts/python/openbabel/openbabel.py
      COMMAND ${SWIG_EXECUTABLE} -python -c++ -small -O -templatereduce -naturalvar -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/python/openbabel/openbabel-python.cpp ${eigen_define} -outdir ${openbabel_SOURCE_DIR}/scripts/python/openbabel ${openbabel_SOURCE_DIR}/scripts/openbabel-python.i
              MAIN_DEPENDENCY openbabel-python.i
              VERBATIM
              )
    endif(RUN_SWIG)
    
    configure_file(${openbabel_SOURCE_DIR}/scripts/python/openbabel/__init__.py.in
                   ${openbabel_SOURCE_DIR}/scripts/python/openbabel/__init__.py)
    add_library(bindings_python MODULE ${openbabel_SOURCE_DIR}/scripts/python/openbabel/openbabel-python.cpp)

    if(APPLE)
      # Don't link against Python library on Mac
      set_target_properties(bindings_python PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
      if(BINDINGS_ONLY)
        target_link_libraries(bindings_python ${BABEL_SYSTEM_LIBRARY})
      else()
        target_link_libraries(bindings_python ${BABEL_LIBRARY})
      endif()
    else()
      # Link against Python library on Linux and Windows
      if(BINDINGS_ONLY)
        target_link_libraries(bindings_python ${PYTHON_LIBRARIES} ${BABEL_SYSTEM_LIBRARY})
      else()
        target_link_libraries(bindings_python ${PYTHON_LIBRARIES} ${BABEL_LIBRARY})
      endif()
    endif()

    if(NOT WIN32)
        set_target_properties(bindings_python PROPERTIES
            OUTPUT_NAME _openbabel
            PREFIX ""
            SUFFIX .so )
        execute_process(
          COMMAND
          ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"
          OUTPUT_VARIABLE PYTHON_INSTDIR
          OUTPUT_STRIP_TRAILING_WHITESPACE
        )
	#workaround for https://bugs.launchpad.net/ubuntu/+source/python3-defaults/+bug/1814653
	if(NOT ${PYTHON_INSTDIR} MATCHES "python[0-9].[0-9]")
	  execute_process(
           COMMAND
           ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(1,1,prefix='${CMAKE_INSTALL_PREFIX}'))"
           OUTPUT_VARIABLE PYTHON_INSTDIR
           OUTPUT_STRIP_TRAILING_WHITESPACE
          )
          set(PYTHON_INSTDIR "${PYTHON_INSTDIR}/dist-packages")
	endif()
        if(NOT BINDINGS_ONLY)
            add_dependencies(bindings_python openbabel)
        endif()

        install(TARGETS bindings_python
                LIBRARY DESTINATION ${PYTHON_INSTDIR}/openbabel
                COMPONENT bindings_python)
        install(FILES ${openbabel_SOURCE_DIR}/scripts/python/openbabel/__init__.py
                DESTINATION ${PYTHON_INSTDIR}/openbabel
                COMPONENT bindings_python)
        install(FILES ${openbabel_SOURCE_DIR}/scripts/python/openbabel/openbabel.py
                DESTINATION ${PYTHON_INSTDIR}/openbabel
                COMPONENT bindings_python)
        install(FILES ${openbabel_SOURCE_DIR}/scripts/python/openbabel/pybel.py
                DESTINATION ${PYTHON_INSTDIR}/openbabel
                COMPONENT bindings_python)
    else(NOT WIN32)
        set_target_properties(bindings_python PROPERTIES
            OUTPUT_NAME _openbabel
            PREFIX ""
            SUFFIX .pyd )
        execute_process(
          COMMAND
          ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"
          OUTPUT_VARIABLE PYTHON_INSTDIR
          OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(NOT BINDINGS_ONLY)
            add_dependencies(bindings_python openbabel)
        endif()

        install(TARGETS bindings_python
                LIBRARY DESTINATION ${PYTHON_INSTDIR}/openbabel
                COMPONENT bindings_python)
        install(FILES ${openbabel_SOURCE_DIR}/scripts/python/openbabel/__init__.py
                DESTINATION ${PYTHON_INSTDIR}/openbabel
                COMPONENT bindings_python)
        install(FILES ${openbabel_SOURCE_DIR}/scripts/python/openbabel/openbabel.py
                DESTINATION ${PYTHON_INSTDIR}/openbabel
                COMPONENT bindings_python)
        install(FILES ${openbabel_SOURCE_DIR}/scripts/python/openbabel/pybel.py
                DESTINATION ${PYTHON_INSTDIR}/openbabel
                COMPONENT bindings_python)
        add_custom_command(TARGET bindings_python POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy ${openbabel_SOURCE_DIR}/scripts/python/openbabel/openbabel.py ${openbabel_BINARY_DIR}/bin/${CMAKE_CFG_INTDIR}
          )
    endif(NOT WIN32)
endif (DO_PYTHON_BINDINGS)

###################
# Java
###################

if (JAVA_BINDINGS)
  if (EXISTS ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp OR RUN_SWIG)

    find_package(JNI)
    if (NOT JNI_INCLUDE_DIRS)
      message(STATUS "Java Native Interface NOT found. Java bindings will NOT be compiled.")
    endif ()

    if (RUN_SWIG)
      find_package(Java)
      if (NOT JAVA_COMPILE)
        message(STATUS "Java compiler NOT found. Java bindings will NOT be generated.")
        set(JAVAC_FAIL ON BOOL)
      endif (NOT JAVA_COMPILE)
    endif (RUN_SWIG)

    if(JNI_INCLUDE_DIRS AND NOT JAVAC_FAIL)
      set(DO_JAVA_BINDINGS ON BOOL)
      message(STATUS "Java bindings will be compiled")
    endif()

  else (EXISTS ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp OR RUN_SWIG)
    message(STATUS "Warning: Java bindings NOT found. Generate using -DRUN_SWIG=ON.")

  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp OR RUN_SWIG)
endif (JAVA_BINDINGS)

if (DO_JAVA_BINDINGS)

include_directories(
	${JNI_INCLUDE_DIRS}
    )

    if (RUN_SWIG)
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp
          COMMAND ${CMAKE_COMMAND} -E make_directory ${openbabel_BINARY_DIR}/scripts/java/org/openbabel
          COMMAND ${SWIG_EXECUTABLE} -java -package org.openbabel -c++ -small -O -templatereduce -naturalvar -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp ${eigen_define} -outdir ${openbabel_BINARY_DIR}/scripts/java/org/openbabel ${openbabel_SOURCE_DIR}/scripts/openbabel-java.i
			 COMMAND ${JAVA_COMPILE} ARGS "java/org/openbabel/*.java"
          COMMAND ${CMAKE_COMMAND} -E chdir ${openbabel_BINARY_DIR}/scripts/java ${JAVA_ARCHIVE} cf openbabel.jar org
          COMMAND ${CMAKE_COMMAND} -E copy ${openbabel_BINARY_DIR}/scripts/java/openbabel.jar ${openbabel_SOURCE_DIR}/scripts/java
          COMMAND ${CMAKE_COMMAND} -E remove_directory ${openbabel_BINARY_DIR}/scripts/java/org
          MAIN_DEPENDENCY openbabel-java.i
      )
    endif (RUN_SWIG)

    add_library(bindings_java MODULE ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp)
    if(BINDINGS_ONLY)
      target_link_libraries(bindings_java ${JNI_LIBRARIES} ${BABEL_SYSTEM_LIBRARY})
    else()
      target_link_libraries(bindings_java ${JNI_LIBRARIES} ${BABEL_LIBRARY})
    endif()

    set_target_properties(bindings_java PROPERTIES
        OUTPUT_NAME openbabel_java)
    if(APPLE)
      set_target_properties(bindings_java PROPERTIES SUFFIX ".jnilib")
    endif(APPLE)
    if(NOT BINDINGS_ONLY)
      add_dependencies(bindings_java openbabel)
    endif()

    install(TARGETS bindings_java
            LIBRARY DESTINATION ${LIB_INSTALL_DIR}
            COMPONENT bindings_java)
    install(FILES ${openbabel_SOURCE_DIR}/scripts/java/openbabel.jar
            DESTINATION ${LIB_INSTALL_DIR}
            COMPONENT bindings_java)

endif (DO_JAVA_BINDINGS)


###################
# R
###################

if (R_BINDINGS)
	if (EXISTS ${openbabel_SOURCE_DIR}/scripts/R/openbabel-R.cpp OR RUN_SWIG)

    if (RUN_SWIG)
      set(DO_R_BINDINGS ON BOOL)
    endif (RUN_SWIG)


  else (EXISTS ${openbabel_SOURCE_DIR}/scripts/R/openbabel-R.cpp OR RUN_SWIG)
	  message(STATUS "Warning: R bindings NOT found. Generate using -DRUN_SWIG=ON.")

  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/R/openbabel-R.cpp OR RUN_SWIG)
endif (R_BINDINGS)

if (DO_R_BINDINGS)


    include_directories(/usr/share/R/include )
    if (RUN_SWIG)
		 add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/R/openbabel-R.cpp
			 COMMAND ${CMAKE_COMMAND} -E make_directory ${openbabel_BINARY_DIR}/scripts/R
			 COMMAND ${SWIG_EXECUTABLE} -r -package openbabelR  -c++ -small -O -templatereduce -naturalvar -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/R/openbabel-R.cpp ${eigen_define} -outdir ${openbabel_BINARY_DIR}/scripts/R ${openbabel_SOURCE_DIR}/scripts/openbabel-R.i
			 # COMMAND sed  -i  -e "'s/reg\.finalizer(\(.*\)) /reg.finalizer(\1); /g'" ${openbabel_BINARY_DIR}/scripts/R/openbabelR.R
			 MAIN_DEPENDENCY openbabel-R.i
      )
    endif (RUN_SWIG)

	 add_library(bindings_R MODULE ${openbabel_SOURCE_DIR}/scripts/R/openbabel-R.cpp)
	 if(BINDINGS_ONLY)
		 target_link_libraries(bindings_R ${BABEL_SYSTEM_LIBRARY})
	 else()
		 target_link_libraries(bindings_R ${BABEL_LIBRARY})
	 endif()

	 set_target_properties(bindings_R PROPERTIES
		 OUTPUT_NAME openbabel_R)
    if(NOT BINDINGS_ONLY)
		 add_dependencies(bindings_R openbabel)
    endif()

	 install(TARGETS bindings_R
            LIBRARY DESTINATION ${LIB_INSTALL_DIR}
				COMPONENT bindings_R)
	 install(FILES ${openbabel_SOURCE_DIR}/scripts/R/openbabel_java.R
            DESTINATION ${LIB_INSTALL_DIR}
				COMPONENT bindings_R)

endif (DO_R_BINDINGS)

###################
# Perl
###################

if (PERL_BINDINGS)
  if (EXISTS ${openbabel_SOURCE_DIR}/scripts/perl/openbabel-perl.cpp OR RUN_SWIG)

    find_package(PerlLibs)
    if (NOT PERL_LIBRARY)
      message(STATUS "Perl library files NOT found. Perl bindings will NOT be compiled.")
    endif (NOT PERL_LIBRARY)

    if(PERL_LIBRARY)
      set(DO_PERL_BINDINGS ON BOOL)
      message(STATUS "Perl bindings will be compiled")
      ### On some Linux OSes, linking against Perl fails due to no -lpthread.
      ### The following is an attempt to work around this by taking the
      ### link libraries from Perl.
      execute_process(
        COMMAND
          ${PERL_EXECUTABLE} -V:perllibs
        OUTPUT_VARIABLE
          PERLLIBS_OUTPUT_VARIABLE
        RESULT_VARIABLE
          PERLLIBS_RESULT_VARIABLE
        )
      if (NOT PERLLIBS_RESULT_VARIABLE)
        string(REGEX REPLACE "perllibs='([^']+)'.*" "\\1" PERL_EXTRA_LINK_FLAGS ${PERLLIBS_OUTPUT_VARIABLE})
      endif (NOT PERLLIBS_RESULT_VARIABLE)
    endif(PERL_LIBRARY)

  else (EXISTS ${openbabel_SOURCE_DIR}/scripts/perl/openbabel-perl.cpp OR RUN_SWIG)
    message(STATUS "Warning: Perl bindings NOT found. Generate using -DRUN_SWIG=ON.")

  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/perl/openbabel-perl.cpp OR RUN_SWIG)
endif (PERL_BINDINGS)

if (DO_PERL_BINDINGS)
    include_directories(
        ${PERL_INCLUDE_PATH}
    )

    if (RUN_SWIG)
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/perl/openbabel-perl.cpp
          COMMAND ${SWIG_EXECUTABLE} -perl5 -c++ -small -O -templatereduce -naturalvar -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/perl/openbabel-perl.cpp -outdir ${openbabel_SOURCE_DIR}/scripts/perl ${openbabel_SOURCE_DIR}/scripts/openbabel-perl.i
          MAIN_DEPENDENCY openbabel-perl.i
      )
      if (NOT WIN32)
        add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/perl/OpenBabel.pm
          COMMAND sed -e \"s/^@EXPORT.*/& sub dl_load_flags { 0x01 }/\" <${openbabel_SOURCE_DIR}/scripts/perl/OpenBabel.pm >${openbabel_SOURCE_DIR}/scripts/perl/OpenBabel.new
          COMMAND ${CMAKE_COMMAND} -E copy ${openbabel_SOURCE_DIR}/scripts/perl/OpenBabel.new ${openbabel_SOURCE_DIR}/scripts/perl/OpenBabel.pm
          COMMAND ${CMAKE_COMMAND} -E remove ${openbabel_SOURCE_DIR}/scripts/perl/OpenBabel.new
          MAIN_DEPENDENCY openbabel-perl.i
        )
      endif (NOT WIN32)
    endif (RUN_SWIG)

    add_library(bindings_perl MODULE ${openbabel_SOURCE_DIR}/scripts/perl/openbabel-perl.cpp)
    if(BINDINGS_ONLY)
      target_link_libraries(bindings_perl ${PERL_LIBRARY} ${BABEL_SYSTEM_LIBRARY} ${PERL_EXTRA_LINK_FLAGS})
    else()
      target_link_libraries(bindings_perl ${PERL_LIBRARY} ${BABEL_LIBRARY} ${PERL_EXTRA_LINK_FLAGS})
    endif()

    set_target_properties(bindings_perl PROPERTIES
        PREFIX ""
        OUTPUT_NAME OpenBabel)
    if(NOT BINDINGS_ONLY)
      add_dependencies(bindings_perl openbabel)
    endif()

    install(TARGETS bindings_perl
            LIBRARY DESTINATION ${LIB_INSTALL_DIR}/auto/Chemistry/OpenBabel/
            COMPONENT bindings_perl)
    install(FILES ${openbabel_SOURCE_DIR}/scripts/perl/OpenBabel.pm
            DESTINATION ${LIB_INSTALL_DIR}/Chemistry/
            COMPONENT bindings_perl)

endif (DO_PERL_BINDINGS)

###################
# Ruby
###################

if (RUBY_BINDINGS)
  if (EXISTS ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.cpp OR RUN_SWIG)
    find_package(Ruby)
    if (NOT RUBY_LIBRARY)
      message(STATUS "Ruby library files NOT found. Ruby bindings will NOT be compiled.")
    endif (NOT RUBY_LIBRARY)

    if (RUBY_LIBRARY)
      set(DO_RUBY_BINDINGS ON BOOL)
      message(STATUS "Ruby bindings will be compiled")
    endif (RUBY_LIBRARY)

  else (EXISTS ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.cpp OR RUN_SWIG)
    message(STATUS "Warning: Ruby bindings NOT found. Generate using -DRUN_SWIG=ON.")

  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.cpp OR RUN_SWIG)
endif (RUBY_BINDINGS)

if (DO_RUBY_BINDINGS)

    if (RUN_SWIG)
      set(WORKAROUND_CMAKE260_BUG "s/rb_define_module(\"Openbabel\")/rb_define_module(\"OpenBabel\")/")
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.cpp
          COMMAND ${SWIG_EXECUTABLE} -ruby -c++ -small -O -templatereduce -naturalvar -autorename -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.txt ${openbabel_SOURCE_DIR}/scripts/openbabel-ruby.i
          COMMAND sed -e '${WORKAROUND_CMAKE_260_BUG}' -e's/void Init_OpenBabel/void Init_openbabel/' <${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.txt >${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.cpp
          COMMAND ${CMAKE_COMMAND} -E remove ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.txt
          MAIN_DEPENDENCY openbabel-ruby.i
      )
    endif (RUN_SWIG)

    include_directories(
        ${RUBY_INCLUDE_PATH}
    )

    add_library(bindings_ruby MODULE ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.cpp)
    if(BINDINGS_ONLY)
      target_link_libraries(bindings_ruby ${RUBY_LIBRARY} ${BABEL_SYSTEM_LIBRARY})
    else()
      target_link_libraries(bindings_ruby ${RUBY_LIBRARY} ${BABEL_LIBRARY})
    endif()

    set_target_properties(bindings_ruby PROPERTIES
        OUTPUT_NAME openbabel
        PREFIX "" )

    if(NOT BINDINGS_ONLY)
      add_dependencies(bindings_ruby openbabel)
    endif()

    install(TARGETS bindings_ruby
            LIBRARY DESTINATION ${LIB_INSTALL_DIR}
            COMPONENT bindings_ruby)

endif (DO_RUBY_BINDINGS)

###################
# CSharp
###################

if (CSHARP_BINDINGS)
  if (EXISTS ${openbabel_SOURCE_DIR}/scripts/csharp/openbabel-csharp.cpp OR RUN_SWIG)
    if (NOT CSHARP_EXECUTABLE)
      message(STATUS "CSHARP_EXECUTABLE NOT defined. .NET bindings will NOT be created.")
    endif (NOT CSHARP_EXECUTABLE)

    if (CSHARP_EXECUTABLE)
      set(DO_CSHARP_BINDINGS ON BOOL)
      message(STATUS "CSharp bindings will be compiled")
    endif (CSHARP_EXECUTABLE)

  else (EXISTS ${openbabel_SOURCE_DIR}/scripts/csharp/openbabel-csharp.cpp OR RUN_SWIG)
    message(STATUS "Warning: Csharp bindings NOT found. Generate using -DRUN_SWIG=ON.")
  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/csharp/openbabel-csharp.cpp OR RUN_SWIG)
endif (CSHARP_BINDINGS)

if (DO_CSHARP_BINDINGS)

    if (RUN_SWIG)
      set(PLATFORM_TYPE "")
      if (WIN32) # i.e. any Windows
          if (CMAKE_GENERATOR MATCHES Win64$)
              set(PLATFORM_TYPE "/platform:anycpu")
          else()
              set(PLATFORM_TYPE "/platform:x86")
          endif()
      endif(WIN32)
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/csharp/openbabel-csharp.cpp ${openbabel_SOURCE_DIR}/scripts/csharp/OBDotNet.dll
          COMMAND ${CMAKE_COMMAND} -E remove_directory ${openbabel_SOURCE_DIR}/scripts/csharp/src
          COMMAND ${CMAKE_COMMAND} -E make_directory ${openbabel_SOURCE_DIR}/scripts/csharp/src
          COMMAND ${CMAKE_COMMAND} -E copy ${openbabel_SOURCE_DIR}/scripts/csharp/OBDotNetAssemblyInfo.cs ${openbabel_SOURCE_DIR}/scripts/csharp/src
          COMMAND ${SWIG_EXECUTABLE} -csharp -c++ -small -O -templatereduce -namespace OpenBabel -outdir ${openbabel_SOURCE_DIR}/scripts/csharp/src -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/csharp/openbabel-csharp.cpp ${openbabel_SOURCE_DIR}/scripts/openbabel-csharp.i
          COMMAND ${CMAKE_COMMAND} -E chdir ${openbabel_SOURCE_DIR}/scripts/csharp/src ${CSHARP_EXECUTABLE} /target:library ${PLATFORM_TYPE} /keyfile:${openbabel_SOURCE_DIR}/scripts/csharp/obdotnet.snk /optimize /out:${openbabel_SOURCE_DIR}/scripts/csharp/OBDotNet.dll *.cs
          COMMAND ${CMAKE_COMMAND} -E remove_directory ${openbabel_SOURCE_DIR}/scripts/csharp/src
          MAIN_DEPENDENCY openbabel-csharp.i
          WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts
      )
    endif (RUN_SWIG)

    add_library(bindings_csharp MODULE ${openbabel_SOURCE_DIR}/scripts/csharp/openbabel-csharp.cpp)
    target_link_libraries(bindings_csharp ${BABEL_LIBRARY})
    set_target_properties(bindings_csharp PROPERTIES
        OUTPUT_NAME openbabel_csharp )
    install(TARGETS bindings_csharp LIBRARY DESTINATION ${LIB_INSTALL_DIR})
    install(FILES ${openbabel_SOURCE_DIR}/scripts/csharp/OBDotNet.dll DESTINATION ${LIB_INSTALL_DIR})

endif (DO_CSHARP_BINDINGS)

###################
# PHP
###################

if (PHP_BINDINGS)
  if (EXISTS ${openbabel_SOURCE_DIR}/scripts/php/openbabel-php.cpp OR RUN_SWIG)

    FIND_PROGRAM(PHP_CONFIG_EXECUTABLE
      NAMES php5-config php-config
    )
    if (NOT PHP_CONFIG_EXECUTABLE)
      message(STATUS "PHP-Config files NOT found. PHP bindings will NOT be compiled.")
    endif (NOT PHP_CONFIG_EXECUTABLE)

    if(PHP_CONFIG_EXECUTABLE)
      set(DO_PHP_BINDINGS ON BOOL)
      message(STATUS "PHP bindings will be compiled")
      execute_process(
        COMMAND
          ${PHP_CONFIG_EXECUTABLE} --includes
        OUTPUT_VARIABLE
          PHP_CONFIG_INCLUDES
        RESULT_VARIABLE
          PHP_CONFIG_RESULT
      )
      STRING(REPLACE "\n" "" PHP_CONFIG_INCLUDES "${PHP_CONFIG_INCLUDES}")
      STRING(REPLACE "-I" "" PHP_CONFIG_INCLUDES "${PHP_CONFIG_INCLUDES}")
      STRING(REPLACE " " ";" PHP_CONFIG_INCLUDES "${PHP_CONFIG_INCLUDES}")

      execute_process(
      	COMMAND
      	  ${PHP_CONFIG_EXECUTABLE} --extension-dir
      	OUTPUT_VARIABLE
      	  PHP_EXTENSION_DIR
      )
      STRING(REPLACE "\n" "" PHP_EXTENSION_DIR "${PHP_EXTENSION_DIR}")
    endif(PHP_CONFIG_EXECUTABLE)

  else (EXISTS ${openbabel_SOURCE_DIR}/scripts/php/openbabel-php.cpp OR RUN_SWIG)
    message(STATUS "Warning: PHP bindings NOT found. Generate using -DRUN_SWIG=ON.")
  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/php/openbabel-php.cpp OR RUN_SWIG)
endif (PHP_BINDINGS)

if (DO_PHP_BINDINGS)
    include_directories(
        ${PHP_CONFIG_INCLUDES}
    )
    if (RUN_SWIG)
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/php/openbabel-php.cpp
          COMMAND ${SWIG_EXECUTABLE} -php -c++ -small -O -templatereduce -naturalvar -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/php/openbabel-php.cpp -outdir ${openbabel_SOURCE_DIR}/scripts/php ${openbabel_SOURCE_DIR}/scripts/openbabel-php.i
          COMMAND sed -i -e's/abstract class OBForceField/class OBForceField/g' -e's/abstract class OBFingerprint/class OBFingerprint/g' -e's/abstract class OBOp/class OBOp/g' -e 's/static function FindType/function FindType/g' ${openbabel_SOURCE_DIR}/scripts/php/openbabel.php
          MAIN_DEPENDENCY openbabel-php.i
      )
    endif (RUN_SWIG)

    # set openbabel.php include path
    execute_process(COMMAND sed -i -e "s,include_once '\\(.*\\)openbabel.php,include_once '${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}/openbabel.php,g" ${openbabel_SOURCE_DIR}/scripts/php/baphpel.php)

    add_library(bindings_php MODULE ${openbabel_SOURCE_DIR}/scripts/php/openbabel-php.cpp)
    # Overriding --no-undefined linker flag, since PHP library should not be linked into extension.
    set(PHP_EXTRA_LINK_FLAGS "-Wl,--unresolved-symbols=ignore-all")
    if(BINDINGS_ONLY)
      target_link_libraries(bindings_php ${BABEL_SYSTEM_LIBRARY} ${PHP_EXTRA_LINK_FLAGS})
    else()
      target_link_libraries(bindings_php ${BABEL_LIBRARY} ${PHP_EXTRA_LINK_FLAGS})
    endif()

    set_target_properties(bindings_php PROPERTIES
        PREFIX ""
        OUTPUT_NAME openbabel-php)
    if(NOT BINDINGS_ONLY)
      add_dependencies(bindings_php openbabel)
    endif()

    install(TARGETS bindings_php
            LIBRARY DESTINATION ${PHP_EXTENSION_DIR}
	    COMPONENT bindings_php)
    install(FILES ${openbabel_SOURCE_DIR}/scripts/php/openbabel.php
            DESTINATION ${LIB_INSTALL_DIR}
	    COMPONENT bindings_php)
	install(FILES ${openbabel_SOURCE_DIR}/scripts/php/baphpel.php
	        DESTINATION ${LIB_INSTALL_DIR}
	    COMPONENT bindings_php)

endif (DO_PHP_BINDINGS)
