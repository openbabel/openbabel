set(SRCS
     OBGUI.cpp stdwx.cpp optswx.cpp)
include_directories(.)
add_executable(OBGUI ${SRCS} ${HDRS})
target_link_libraries(OBGUI ${BABEL_LIBRARY} ${wxWidgets_LIBRARIES})
install(TARGETS OBGUI 
	  RUNTIME DESTINATION bin
	  LIBRARY DESTINATION lib
	  ARCHIVE DESTINATION lib
)
