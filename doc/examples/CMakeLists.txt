set(examples
  obconformersearch_default
  obconversion_readstring
)

foreach(example ${examples})
  add_executable(${example} ${example}.cpp)
  target_link_libraries(${example} ${BABEL_LIBRARY})
  install(TARGETS ${example}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )
endforeach(example)
