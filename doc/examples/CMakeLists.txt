set(examples
    obconformersearch_default
    obconversion_readstring

)

foreach(example ${examples})
  add_executable(${example} ${example}.cpp)
  target_link_libraries(${example} ${BABEL_LIBRARY})
  if(MSVC)
    set_target_properties(${example} PROPERTIES
	    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endif(MSVC)
  install(TARGETS ${example}
                  RUNTIME DESTINATION bin
                  LIBRARY DESTINATION lib
                  ARCHIVE DESTINATION lib
  )
endforeach(example)
