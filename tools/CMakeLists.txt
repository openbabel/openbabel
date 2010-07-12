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
endif(MINIMAL_BUILD)

foreach(tool ${tools})
  add_executable(${tool} ${tool}.cpp)
  target_link_libraries(${tool} ${BABEL_LIBRARY})
  if(MSVC)
    set_target_properties(${tool} PROPERTIES
	    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endif(MSVC)
  install(TARGETS ${tool}
                  RUNTIME DESTINATION bin
                  LIBRARY DESTINATION lib
                  ARCHIVE DESTINATION lib
  )
endforeach(tool)

if(NOT MINIMAL_BUILD)
  # obspectrophore -- requires getopt
  add_executable(obspectrophore obspectrophore.cpp getopt.c)
  target_link_libraries(obspectrophore ${BABEL_LIBRARY})
  if(MSVC)
    set_target_properties(obspectrophore PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endif(MSVC)
  install(TARGETS obspectrophore
                  RUNTIME DESTINATION bin
                  LIBRARY DESTINATION lib
                  ARCHIVE DESTINATION lib
  )

  # obgrep -- requires getopt
  add_executable(obgrep obgrep.cpp getopt.c)
  target_link_libraries(obgrep ${BABEL_LIBRARY})
  if(MSVC)
    set_target_properties(obgrep PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endif(MSVC)
  install(TARGETS obgrep
                  RUNTIME DESTINATION bin
                  LIBRARY DESTINATION lib
                  ARCHIVE DESTINATION lib
  )

endif(NOT MINIMAL_BUILD)
