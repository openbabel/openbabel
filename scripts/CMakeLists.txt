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

if (EIGEN2_FOUND)
  set (eigen2_define "-DHAVE_EIGEN2")
  set (CPPFLAGS_VALUE "-I${EIGEN2_INCLUDE_DIR}")
  set (CPPFLAGS_TEXT "CPPFLAGS=${CPPFLAGS_VALUE}")
endif (EIGEN2_FOUND)

###################
# Python
###################

if (PYTHON_BINDINGS)
  if (EXISTS ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp OR RUN_SWIG)
    find_package(PythonLibs)
    if (NOT PYTHONLIBS_FOUND)
      message(STATUS "Python libraries NOT found")
    endif (NOT PYTHONLIBS_FOUND)

    find_package(PythonInterp)
    if (NOT PYTHONINTERP_FOUND)
      message(STATUS "Python interpreter NOT found")
    endif (NOT PYTHONINTERP_FOUND)

    if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
      set(DO_PYTHON_BINDINGS ON BOOL)
      message(STATUS "Python bindings will be compiled")
    endif(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)

  else (EXISTS ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp OR RUN_SWIG)
    message(STATUS "Warning: Python bindings NOT found. Generate using -DRUN_SWIG.")

  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp OR RUN_SWIG)
endif (PYTHON_BINDINGS)

if (DO_PYTHON_BINDINGS)

    include_directories(
        ${PYTHON_INCLUDE_PATH}
    )

    if (RUN_SWIG)
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp ${openbabel_SOURCE_DIR}/scripts/python/openbabel.py
          COMMAND ${SWIG_EXECUTABLE} -python -c++ -small -O -templatereduce -naturalvar -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp ${eigen2_define} -outdir ${openbabel_SOURCE_DIR}/scripts/python ${openbabel_SOURCE_DIR}/scripts/openbabel-python.i
          MAIN_DEPENDENCY openbabel-python.i
      )
    endif (RUN_SWIG)

    if(NOT WIN32)
        add_custom_target(_openbabel ALL
            COMMAND ${CPPFLAGS_TEXT} ${PYTHON_EXECUTABLE} ${openbabel_SOURCE_DIR}/scripts/python/setup.py build --build-lib=${openbabel_BINARY_DIR}/scripts/pybuild
            DEPENDS openbabel ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp
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
        if (APPLE)
           # Update the built version
          add_custom_command(TARGET _openbabel POST_BUILD
           COMMAND install_name_tool -change ${openbabel_BINARY_DIR}/lib/libopenbabel.${SOVERSION}.dylib ${CMAKE_INSTALL_PREFIX}/lib/libopenbabel.${SOVERSION}.dylib ${openbabel_BINARY_DIR}/scripts/pybuild/_openbabel.so
           VERBATIM)
        endif()
        add_dependencies(_openbabel openbabel)

        if (PYTHON_PREFIX)
            # In an execute_process, we cannot set the CPPFLAGS variable as we did with COMMAND above.
            # ...hence the different approach here.
            install(CODE "set(ENV{CPPFLAGS} ${CPPFLAGS_VALUE})
                          execute_process(COMMAND ${PYTHON_EXECUTABLE} ${openbabel_SOURCE_DIR}/scripts/python/setup.py install --root=\$ENV{DESTDIR}/ --prefix=${PYTHON_PREFIX} WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts)"   LIBRARY DESTINATION ${LIB_INSTALL_DIR})
        else (PYTHON_PREFIX)
            if (${CMAKE_INSTALL_PREFIX} STREQUAL "/usr/local")
                install(CODE "set(ENV{CPPFLAGS} ${CPPFLAGS_VALUE})
                              execute_process(COMMAND ${PYTHON_EXECUTABLE} ${openbabel_SOURCE_DIR}/scripts/python/setup.py install --root=\$ENV{DESTDIR}/ WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts)"   LIBRARY DESTINATION ${LIB_INSTALL_DIR})
           else (${CMAKE_INSTALL_PREFIX} STREQUAL "/usr/local")
                install(CODE "set(ENV{CPPFLAGS} ${CPPFLAGS_VALUE})
                              execute_process(COMMAND ${PYTHON_EXECUTABLE} ${openbabel_SOURCE_DIR}/scripts/python/setup.py install --root=\$ENV{DESTDIR}/ --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts)"   LIBRARY DESTINATION ${LIB_INSTALL_DIR})
            endif ()
        endif ()
        if (APPLE)
           EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c
               "from sys import stdout; from distutils import sysconfig; stdout.write(sysconfig.get_python_lib(1))"
           OUTPUT_VARIABLE PYTHON_LIB_PATH)
           install(CODE "execute_process(COMMAND install_name_tool -change ${openbabel_BINARY_DIR}/lib/libopenbabel.${SOVERSION}.dylib ${CMAKE_INSTALL_PREFIX}/lib/libopenbabel.${SOVERSION}.dylib ${PYTHON_LIB_PATH}/_openbabel.so)")
        endif()

    else(NOT WIN32)
        add_library(bindings_python MODULE ${openbabel_SOURCE_DIR}/scripts/python/openbabel-python.cpp)
        target_link_libraries(bindings_python ${PYTHON_LIBRARIES} ${BABEL_LIBRARY})
        set_target_properties(bindings_python PROPERTIES
            OUTPUT_NAME _openbabel
            PREFIX ""
            SUFFIX .pyd )
        add_custom_command(TARGET bindings_python POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy ${openbabel_SOURCE_DIR}/scripts/python/openbabel.py ${openbabel_BINARY_DIR}/bin/${CMAKE_CFG_INTDIR}
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
    message(STATUS "Warning: Java bindings NOT found. Generate using -DRUN_SWIG.")

  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp OR RUN_SWIG)
endif (JAVA_BINDINGS)

if (DO_JAVA_BINDINGS)

    include_directories(
        ${JNI_INCLUDE_DIRS}
    )

    if (RUN_SWIG)
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp # ${openbabel_SOURCE_DIR}/scripts/python/openbabel.py
          COMMAND ${CMAKE_COMMAND} -E make_directory ${openbabel_SOURCE_DIR}/scripts/java/org/openbabel
          COMMAND ${SWIG_EXECUTABLE} -java -package org.openbabel -c++ -small -O -templatereduce -naturalvar -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp ${eigen2_define} -outdir ${openbabel_SOURCE_DIR}/scripts/java/org/openbabel ${openbabel_SOURCE_DIR}/scripts/openbabel-java.i
          COMMAND ${CMAKE_COMMAND} -E chdir ${openbabel_SOURCE_DIR}/scripts/java ${JAVA_COMPILE} org/openbabel/*.java
          COMMAND ${CMAKE_COMMAND} -E chdir ${openbabel_SOURCE_DIR}/scripts/java ${JAVA_ARCHIVE} cf openbabel.jar org
          COMMAND ${CMAKE_COMMAND} -E remove_directory ${openbabel_SOURCE_DIR}/scripts/java/org
          MAIN_DEPENDENCY openbabel-java.i
      )
    endif (RUN_SWIG)

    add_library(bindings_java MODULE ${openbabel_SOURCE_DIR}/scripts/java/openbabel-java.cpp)
    target_link_libraries(bindings_java ${JNI_LIBRARIES} ${BABEL_LIBRARY})
    set_target_properties(bindings_java PROPERTIES
        OUTPUT_NAME openbabel_java )

endif (DO_JAVA_BINDINGS)

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
    endif(PERL_LIBRARY)

  else (EXISTS ${openbabel_SOURCE_DIR}/scripts/perl/openbabel-perl.cpp OR RUN_SWIG)
    message(STATUS "Warning: Perl bindings NOT found. Generate using -DRUN_SWIG.")

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
          MAIN_DEPENDENCY ${openbabel_SOURCE_DIR}/scripts/perl/openbabel-perl.cpp
        )
      endif (NOT WIN32)
    endif (RUN_SWIG)

    set(PERL_PREFIX_TEXT "")
    if (OBPERL_PREFIX) # The PERL_PREFIX variable is already taken (it's set by FindPerlLibs :-/ )
      set(PERL_PREFIX_TEXT "PREFIX=${OBPERL_PREFIX}")
    else (OBPERL_PREFIX)
      if (NOT ${CMAKE_INSTALL_PREFIX} STREQUAL "/usr/local")
        set(PERL_PREFIX_TEXT "PREFIX=${CMAKE_INSTALL_PREFIX}")
      endif (NOT ${CMAKE_INSTALL_PREFIX} STREQUAL "/usr/local")
    endif (OBPERL_PREFIX)

    set(PERL_INSTALLDIRS_TEXT "")
    if (OBPERL_INSTALLDIRS)
      set(PERL_INSTALLDIRS_TEXT "INSTALLDIRS=${OBPERL_INSTALLDIRS}")
    endif (OBPERL_INSTALLDIRS)

    if (NOT MSVC)
      add_custom_target(pre_bindings_perl ALL
        # Create the scripts/perl directory; we cannot do this inside
        # bindings_perl as the WORKING_DIRECTORY needs to already exist!
        COMMAND ${CMAKE_COMMAND} -E make_directory ${openbabel_BINARY_DIR}/scripts/perl
        )
      add_custom_target(bindings_perl ALL
        COMMAND ${CMAKE_COMMAND} -E copy ${openbabel_SOURCE_DIR}/scripts/perl/Makefile.PL ${openbabel_BINARY_DIR}/scripts/perl
        COMMAND ${CMAKE_COMMAND} -E copy ${openbabel_SOURCE_DIR}/scripts/perl/OpenBabel.pm ${openbabel_BINARY_DIR}/scripts/perl
        COMMAND ${CMAKE_COMMAND} -E copy ${openbabel_SOURCE_DIR}/scripts/perl/openbabel-perl.cpp ${openbabel_BINARY_DIR}/scripts/perl
        COMMAND PERL5LIB=${openbabel_SOURCE_DIR}/scripts/perl/inc SRC_DIR=${openbabel_SOURCE_DIR} OBJ_DIR=${openbabel_BINARY_DIR} ${PERL_EXECUTABLE} ${openbabel_BINARY_DIR}/scripts/perl/Makefile.PL ${PERL_PREFIX_TEXT} ${PERL_INSTALLDIRS_TEXT}
        COMMAND PERL5LIB=${openbabel_SOURCE_DIR}/scripts/perl/inc make
            DEPENDS openbabel ${openbabel_SOURCE_DIR}/scripts/perl/openbabel-perl.cpp
            WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts/perl
        )
      add_dependencies(bindings_perl pre_bindings_perl openbabel)
      install(CODE "execute_process(COMMAND make install DESTDIR=\$ENV{DESTDIR} WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts/perl)")

    else (NOT MSVC)
      add_library(bindings_perl MODULE ${openbabel_SOURCE_DIR}/scripts/perl/openbabel-perl.cpp)
      target_link_libraries(bindings_perl ${PERL_LIBRARY} ${BABEL_LIBRARY})
      set_target_properties(bindings_perl PROPERTIES
          OUTPUT_NAME openbabel_perl )
    endif(NOT MSVC)

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
    message(STATUS "Warning: Ruby bindings NOT found. Generate using -DRUN_SWIG.")

  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.cpp OR RUN_SWIG)
endif (RUBY_BINDINGS)

if (DO_RUBY_BINDINGS)

    if (RUN_SWIG)
      set(WORKAROUND_CMAKE260_BUG "s/rb_define_module(\"Openbabel\")/rb_define_module(\"OpenBabel\")/")
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.cpp
          COMMAND ${SWIG_EXECUTABLE} -ruby -c++ -small -O -templatereduce -naturalvar -autorename -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.txt ${openbabel_SOURCE_DIR}/scripts/openbabel-ruby.i
          COMMAND sed -e '${WORKAROUND_CMAKE_260_BUG}' -e's/void Init_OpenBabel/void Init_openbabel/' <${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.txt >${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.cpp 
          MAIN_DEPENDENCY openbabel-ruby.i
      )
    endif (RUN_SWIG)

    include_directories(
        ${RUBY_INCLUDE_PATH}
    )

    if(NOT MSVC)
      add_custom_target(pre_bindings_ruby ALL
        # Create the scripts/ruby directory; we cannot do this inside
        # bindings_ruby as the WORKING_DIRECTORY needs to already exist!
        COMMAND ${CMAKE_COMMAND} -E make_directory ${openbabel_BINARY_DIR}/scripts/ruby
        )
      set (DESTDIR_TEXT "\$ENV{DESTDIR}")
      if (RUBY_PREFIX)
        # It is reasonably safe to assume that RUBY_POSSIBLE_LIB_DIR contains
        # no regex metacharacters except possibly '.':
        string(REGEX REPLACE ^${RUBY_POSSIBLE_LIB_DIR}/ruby ${RUBY_PREFIX} ruby_install_dir "${RUBY_SITEARCH_DIR}")
      else (RUBY_PREFIX)
        execute_process(COMMAND ${RUBY_EXECUTABLE} -r rbconfig -e "print Config::CONFIG['exec_prefix']"
          OUTPUT_VARIABLE ruby_exec_prefix)
        # It is reasonably safe to assume that ruby_exec_prefix contains no
        # regex metacharacters except possibly '.':
        string(REGEX REPLACE ^${ruby_exec_prefix} ${CMAKE_INSTALL_PREFIX} ruby_install_dir "${RUBY_SITEARCH_DIR}")
      endif (RUBY_PREFIX)
      add_custom_target(bindings_ruby ALL
        COMMAND ${RUBY_EXECUTABLE} ${openbabel_SOURCE_DIR}/scripts/ruby/extconf.rb --with-openbabel-lib=${openbabel_BINARY_DIR}/lib${LIB_SUFFIX} --with-openbabel-include=${openbabel_BINARY_DIR}/include
        COMMAND ${CMAKE_COMMAND} -E echo "Compiling Ruby bindings"
        COMMAND make
        DEPENDS ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.cpp
        WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts/ruby
        )
      add_dependencies(bindings_ruby pre_bindings_ruby openbabel)
      install(CODE "execute_process(COMMAND make install DESTDIR=\$ENV{DESTDIR} WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts/ruby)")
   else(NOT MSVC)
      add_library(bindings_ruby MODULE ${openbabel_SOURCE_DIR}/scripts/ruby/openbabel-ruby.cpp)
      target_link_libraries(bindings_ruby ${RUBY_LIBRARY} ${BABEL_LIBRARY})
      set_target_properties(bindings_ruby PROPERTIES
          OUTPUT_NAME openbabel_ruby )
    endif(NOT MSVC)

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
    message(STATUS "Warning: Csharp bindings NOT found. Generate using -DRUN_SWIG.")
  endif (EXISTS ${openbabel_SOURCE_DIR}/scripts/csharp/openbabel-csharp.cpp OR RUN_SWIG)
endif (CSHARP_BINDINGS)

if (DO_CSHARP_BINDINGS)

    if (RUN_SWIG)
      set(PLATFORM_TYPE "")
      if (WIN32) # We need to enforce 32-bit .NET dlls or it won't load
          set(PLATFORM_TYPE "/platform:x86")
      endif (WIN32)
      add_custom_command(OUTPUT ${openbabel_SOURCE_DIR}/scripts/csharp/openbabel-csharp.cpp ${openbabel_SOURCE_DIR}/scripts/csharp/OBDotNet.dll
          COMMAND ${CMAKE_COMMAND} -E remove_directory ${openbabel_BINARY_DIR}/scripts/csharp
          COMMAND ${CMAKE_COMMAND} -E make_directory ${openbabel_BINARY_DIR}/scripts/csharp
          COMMAND ${CMAKE_COMMAND} -E copy ${openbabel_SOURCE_DIR}/windows-vc2008/Distribution/OBDotNetAssemblyInfo.cs ${openbabel_BINARY_DIR}/scripts/csharp
          COMMAND ${SWIG_EXECUTABLE} -csharp -c++ -small -O -templatereduce -namespace OpenBabel -outdir ${openbabel_BINARY_DIR}/scripts/csharp -I${openbabel_SOURCE_DIR}/include -I${openbabel_BINARY_DIR}/include -o ${openbabel_SOURCE_DIR}/scripts/csharp/openbabel-csharp.cpp ${openbabel_SOURCE_DIR}/scripts/openbabel-csharp.i
          COMMAND ${CMAKE_COMMAND} -E chdir ${openbabel_BINARY_DIR}/scripts/csharp ${CSHARP_EXECUTABLE} /target:library ${PLATFORM_TYPE} /keyfile:${openbabel_SOURCE_DIR}/windows-vc2008/Distribution/obdotnet.snk /optimize /out:${openbabel_SOURCE_DIR}/scripts/csharp/OBDotNet.dll *.cs
          MAIN_DEPENDENCY openbabel-csharp.i
          WORKING_DIRECTORY ${openbabel_BINARY_DIR}/scripts
      )
    endif (RUN_SWIG)

    add_library(bindings_csharp MODULE ${openbabel_SOURCE_DIR}/scripts/csharp/openbabel-csharp.cpp)
    target_link_libraries(bindings_csharp ${BABEL_LIBRARY})
    set_target_properties(bindings_csharp PROPERTIES
        OUTPUT_NAME openbabel_csharp )

endif (DO_CSHARP_BINDINGS)

