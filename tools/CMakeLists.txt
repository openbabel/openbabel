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

  # no getopt on windows
  if(NOT WIN32)
    set(tools
        ${tools}
        obgrep
    )
  endif(NOT WIN32)

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
