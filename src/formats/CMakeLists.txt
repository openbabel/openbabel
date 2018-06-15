include(formats.cmake)

if(NOT WIN32)
  set(libs ${libs} m)
endif(NOT WIN32)

if(WIN32)
  set(openbabel_srcs ${openbabel_srcs}
      dlhandler_win32
      )
else(WIN32)
  set(openbabel_srcs ${openbabel_srcs}
      dlhandler_unix
      )
endif(WIN32)


if(MSVC)
  foreach(formatgroup formats_common formats_utility formats_compchem formats_misc ${optional_formatgroups})
    set(formatsources "")
    foreach(format ${${formatgroup}})
      set(formatsources ${formatsources} ${format}.cpp ${${format}_additional_sources})
    endforeach(format)
    add_library(${formatgroup} ${PLUGIN_TYPE} ${formatsources}
        "${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h")
    target_link_libraries(${formatgroup} ${libs} ${XDR_LIBRARY} openbabel)
    install(TARGETS ${formatgroup}
                  RUNTIME DESTINATION ${BIN_INSTALL_DIR}
                  LIBRARY DESTINATION ${OB_PLUGIN_INSTALL_DIR}
                  ARCHIVE DESTINATION ${OB_PLUGIN_INSTALL_DIR})
    set_target_properties(${formatgroup} PROPERTIES
                  OUTPUT_NAME ${formatgroup}
                  PREFIX ""
                  SUFFIX ${MODULE_EXTENSION})
  endforeach(formatgroup)
else(MSVC)
foreach(format ${formats})
  add_library(${format} ${PLUGIN_TYPE} ${format}.cpp
    ${${format}_additional_sources}
    "${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h")
  target_link_libraries(${format} ${libs} openbabel)
  install(TARGETS ${format}
                  RUNTIME DESTINATION ${BIN_INSTALL_DIR}
                  LIBRARY DESTINATION ${OB_PLUGIN_INSTALL_DIR}
                  ARCHIVE DESTINATION ${OB_PLUGIN_INSTALL_DIR})
  set_target_properties(${format} PROPERTIES
                        OUTPUT_NAME ${format}
                        PREFIX ""
                        SUFFIX ${MODULE_EXTENSION})
endforeach(format)
endif(MSVC)

if(NOT MINIMAL_BUILD)
if(MSVC)
  set(xmlsources "")
  foreach(format ${formats_xml})
    set(xmlsources ${xmlsources} xml/${format}.cpp)
  endforeach(format)
  add_library(formats_xml ${PLUGIN_TYPE} ${xmlsources} xml/xml.cpp)
  target_link_libraries(formats_xml ${libs} ${LIBXML2_LIBRARIES} openbabel)
  install(TARGETS formats_xml
                  RUNTIME DESTINATION ${BIN_INSTALL_DIR}
                  LIBRARY DESTINATION ${OB_PLUGIN_INSTALL_DIR}
                  ARCHIVE DESTINATION ${OB_PLUGIN_INSTALL_DIR})
  set_target_properties(formats_xml PROPERTIES
                  OUTPUT_NAME formats_xml
                  PREFIX ""
                  SUFFIX ${MODULE_EXTENSION})

  if(WITH_JSON)
    set(jsonsources "")
    foreach(format ${formats_json})
      set(jsonsources ${jsonsources} json/${format}.cpp)
    endforeach(format)
    add_library(formats_json ${PLUGIN_TYPE} ${jsonsources}
                "${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h")
    target_link_libraries(formats_json ${libs} ${JSON_LIBRARY} openbabel)
    install(TARGETS formats_json
            RUNTIME DESTINATION ${BIN_INSTALL_DIR}
            LIBRARY DESTINATION ${OB_PLUGIN_INSTALL_DIR}
            ARCHIVE DESTINATION ${OB_PLUGIN_INSTALL_DIR})
    set_target_properties(formats_json PROPERTIES
                          OUTPUT_NAME formats_json
                          PREFIX ""
                          SUFFIX ${MODULE_EXTENSION})
  endif()

else(MSVC)
foreach(format ${formats_xml})
  add_library(${format} ${PLUGIN_TYPE} xml/${format}.cpp xml/xml.cpp)
  target_link_libraries(${format} ${libs} ${LIBXML2_LIBRARIES} openbabel)
  install(TARGETS ${format}
                  RUNTIME DESTINATION ${BIN_INSTALL_DIR}
                  LIBRARY DESTINATION ${OB_PLUGIN_INSTALL_DIR}
                  ARCHIVE DESTINATION ${OB_PLUGIN_INSTALL_DIR})
  set_target_properties(${format} PROPERTIES
                        OUTPUT_NAME ${format}
                        PREFIX ""
                        SUFFIX ${MODULE_EXTENSION})
endforeach(format)

if(WITH_JSON)
foreach(format ${formats_json})
  add_library(${format} ${PLUGIN_TYPE} json/${format}.cpp
              "${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h")
  target_link_libraries(${format} ${libs} openbabel)
  install(TARGETS ${format}
          RUNTIME DESTINATION ${BIN_INSTALL_DIR}
          LIBRARY DESTINATION ${OB_PLUGIN_INSTALL_DIR}
          ARCHIVE DESTINATION ${OB_PLUGIN_INSTALL_DIR})
  set_target_properties(${format} PROPERTIES
                        OUTPUT_NAME ${format}
                        PREFIX ""
                        SUFFIX ${MODULE_EXTENSION})
endforeach(format)
endif()
endif(MSVC)
endif(NOT MINIMAL_BUILD)
