if(MINIMAL_BUILD)
  set(tools
        babel
  )
else(MINIMAL_BUILD)
  set(tools
        babel
        obabel
        obchiral
        obconformer
        obenergy
        obfit
        obgen
        obminimize
        obprobe
        obprop
        obrotamer
        obrotate
     )
  if(EIGEN3_FOUND)
    set(tools
        ${tools}
        obrms
       )
  endif()
endif(MINIMAL_BUILD)

INCLUDE(CheckFunctionExists)
SET(GETOPT getopt)
CHECK_FUNCTION_EXISTS(${GETOPT} GETOPT_FOUND)

if(BUILD_SHARED)
  foreach(tool ${tools})
    add_executable(${tool} ${tool}.cpp)
    target_link_libraries(${tool} ${BABEL_LIBRARY})
    install(TARGETS ${tool}
                    RUNTIME DESTINATION bin
                    LIBRARY DESTINATION lib
                    ARCHIVE DESTINATION lib
    )
  endforeach(tool)

  if(NOT MINIMAL_BUILD)
    # obspectrophore -- requires getopt
    if (NOT GETOPT_FOUND)
      add_executable(obspectrophore obspectrophore.cpp getopt.c)
    else(NOT GETOPT_FOUND)
      add_executable(obspectrophore obspectrophore.cpp)
    endif(NOT GETOPT_FOUND)
    target_link_libraries(obspectrophore ${BABEL_LIBRARY})
    install(TARGETS obspectrophore
                    RUNTIME DESTINATION bin
                    LIBRARY DESTINATION lib
                    ARCHIVE DESTINATION lib
    )

    # obgrep -- requires getopt
    if (NOT GETOPT_FOUND)
      add_executable(obgrep obgrep.cpp getopt.c)
    else(NOT GETOPT_FOUND)
      add_executable(obgrep obgrep.cpp)
    endif(NOT GETOPT_FOUND)
      target_link_libraries(obgrep ${BABEL_LIBRARY})
      install(TARGETS obgrep
                      RUNTIME DESTINATION bin
                      LIBRARY DESTINATION lib
                      ARCHIVE DESTINATION lib
      )
  endif(NOT MINIMAL_BUILD)

else(BUILD_SHARED)
###############################################################################
#
# Build STATIC babel with all plugins taken from the static libopenbabel.a
#
###############################################################################

  # Set libs to imported openbabel target
  set(libs openbabel)
  if(WITH_STATIC_INCHI)
    set(libs ${libs} ${LIBXML2_LIBRARIES})
  endif()
  if(WITH_STATIC_LIBXML)
    set(libs ${libs} inchi)
  endif()

  if(NOT MSVC)
    find_package(Threads)
    if(CMAKE_USE_PTHREADS_INIT)
      if(BUILD_MIXED)
        set(libs ${libs} ${CMAKE_THREAD_LIBS_INIT})
      else()
        # Make sure to find the static library with .a suffix (see CMAKE_FIND_LIBRARY_SUFFIXES)
        find_library(PTHREAD_LIBRARY NAMES pthread)
        if(PTHREAD_LIBRARY)
          set(libs ${libs} ${PTHREAD_LIBRARY})
        endif()
      endif()
    endif()

    if(BUILD_MIXED)
      if(CMAKE_COMPILER_IS_GNUCXX)
        if(CMAKE_SYSTEM_NAME MATCHES Linux)

          #
          # Relocatable binaries on linux using -static-libgcc. (GCC only)
          #
          set(CMAKE_CXX_FLAGS "-static-libgcc ${CMAKE_CXX_FLAGS}")
          # make sure the linker will find the static libstdc++
          file(REMOVE "${CMAKE_BINARY_DIR}/libstdc++.a")
          execute_process(COMMAND "g++" "-print-file-name=libstdc++.a"
                          OUTPUT_VARIABLE libstdcpp
                          OUTPUT_STRIP_TRAILING_WHITESPACE)
          execute_process(COMMAND "ln" "-s" "${libstdcpp}"
                          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

          link_directories(${CMAKE_BINARY_DIR})

        endif(CMAKE_SYSTEM_NAME MATCHES Linux)
      endif(CMAKE_COMPILER_IS_GNUCXX)
    else(BUILD_MIXED)
      #
      # Default static build using GNU-like compiler.
      #
      if (NOT APPLE)
            set(CMAKE_CXX_FLAGS "-static ${CMAKE_CXX_FLAGS}")
      endif()
      # APPLE doesn't use -static
    endif(BUILD_MIXED)
  endif()

  # build a static babel only
  add_executable(babel babel.cpp)
  target_link_libraries(babel ${libs})

  if(NOT BUILD_MIXED)
    # Prevent -Wl,-Bdynamic from being added to the end of the link line.
    set_target_properties(babel PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  endif(NOT BUILD_MIXED)

  install(TARGETS babel
          RUNTIME DESTINATION bin
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)

endif(BUILD_SHARED)
