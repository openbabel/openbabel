set(tools
      babel
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

foreach(tool ${tools})
  add_executable(${tool} ${tool}.cpp)
  target_link_libraries(${tool} ${BABEL_LIBRARY})
  install(TARGETS ${tool}
                  RUNTIME DESTINATION bin
                  LIBRARY DESTINATION lib
                  ARCHIVE DESTINATION lib
  )
endforeach(tool)
