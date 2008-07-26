
add_executable(roundtrip roundtrip.cpp)
target_link_libraries(roundtrip ${BABEL_LIBRARY})
install(TARGETS roundtrip
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
)
