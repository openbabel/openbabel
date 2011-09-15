set(SRCS
     OBGUI.cpp stdwx.cpp optswx.cpp OBGUI.rc)
include_directories(.)
add_executable(obgui WIN32 ${SRCS} ${HDRS})
if(WIN32)
  set(wxWidgets_LIBRARIES ${wxWidgets_LIBRARIES} comctl32.lib)
endif(WIN32) 
target_link_libraries(obgui ${BABEL_LIBRARY} ${wxWidgets_LIBRARIES})
install(TARGETS obgui 
	  RUNTIME DESTINATION bin
	  LIBRARY DESTINATION lib
	  ARCHIVE DESTINATION lib
)
