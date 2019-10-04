if(MINIMAL_BUILD)
  set(tools
        obabel
  )
else(MINIMAL_BUILD)
  find_package(Boost COMPONENTS program_options iostreams)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
  endif()


  set(tools
        obabel
        obconformer
        obenergy
        obfit
	obfitall
        obgen
        obminimize
        obmm
        obprobe
        obprop
        obrotamer
        obrotate
        obsym
        obtautomer
        obthermo
     )
  if(EIGEN3_FOUND)
    set(tools
        ${tools}
        obdistgen
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
                    LIBRARY DESTINATION lib${LIB_SUFFIX}
                    ARCHIVE DESTINATION lib${LIB_SUFFIX}
    )
  endforeach(tool)

  if(NOT MINIMAL_BUILD)
    # obgrep, obrms, obspectrophore -- require getopt
    set(toolnames obgrep obspectrophore)
    if(EIGEN3_FOUND)
      set(toolnames ${toolnames} obrms)
    endif()
    foreach(tool ${toolnames})
      if (NOT GETOPT_FOUND)
        add_executable(${tool} ${tool}.cpp getopt.c)
      else(NOT GETOPT_FOUND)
        add_executable(${tool} ${tool}.cpp)
      endif(NOT GETOPT_FOUND)
      target_link_libraries(${tool} ${BABEL_LIBRARY})
      install(TARGETS ${tool}
                      RUNTIME DESTINATION bin
                      LIBRARY DESTINATION lib${LIB_SUFFIX}
                      ARCHIVE DESTINATION lib${LIB_SUFFIX}
      )
    endforeach()

  endif(NOT MINIMAL_BUILD)

else(BUILD_SHARED)
###############################################################################
#
# Build STATIC obabel with all plugins taken from the static libopenbabel.a
#
###############################################################################

  # Set libs to imported openbabel target
  set(libs openbabel)
  if(WITH_STATIC_INCHI)
    set(libs ${libs} ${INCHI_LIBRARY})
  endif()
  if(WITH_STATIC_LIBXML)
    set(libs ${libs} ${LIBXML2_LIBRARIES})
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

  # build a static obabel only
  add_executable(obabel obabel.cpp)
  target_link_libraries(obabel ${libs})

  if(NOT BUILD_MIXED)
    # Prevent -Wl,-Bdynamic from being added to the end of the link line.
    set_target_properties(obabel PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  endif(NOT BUILD_MIXED)

  install(TARGETS obabel
          RUNTIME DESTINATION bin
          LIBRARY DESTINATION lib${LIB_SUFFIX}
          ARCHIVE DESTINATION lib${LIB_SUFFIX})

endif(BUILD_SHARED)
