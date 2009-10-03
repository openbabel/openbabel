set(SRCS
     OBGUI.cpp stdwx.cpp optswx.cpp OBGUI.rc)
include_directories(.)
add_executable(OBGUI WIN32 ${SRCS} ${HDRS})
if(WIN32)
  set(wxWidgets_LIBRARIES ${wxWidgets_LIBRARIES} comctl32.lib)
  set_target_properties(OBGUI PROPERTIES
	  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif(WIN32) 
target_link_libraries(OBGUI ${BABEL_LIBRARY} ${wxWidgets_LIBRARIES})
install(TARGETS OBGUI 
	  RUNTIME DESTINATION bin
	  LIBRARY DESTINATION lib
	  ARCHIVE DESTINATION lib
)
