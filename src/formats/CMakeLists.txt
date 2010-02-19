set(formats_common
  smilesformat
  mdlformat
  mol2format
  pdbformat
  )
set(formats_utility
  copyformat
  MNAformat
  molreport
  nulformat
  pngformat
  povrayformat
  reportformat
  svgformat
  textformat
  titleformat
  )
set(formats_compchem
      adfformat
      dmolformat
      fchkformat
      gamessformat
      gaussformat
      gausscubeformat
      gausszmatformat
      gulpformat
      hinformat
      jaguarformat
      molproformat
      mopacformat
      nwchemformat
      qchemformat
      turbomoleformat
      vaspformat
      zindoformat
      )
set(formats_misc
      APIInterface
      CSRformat
      PQSformat
      MCDLformat
      alchemyformat
      acrformat
      amberformat
      balstformat
      bgfformat
      boxformat
      cacaoformat
      cacheformat
      carformat
      cccformat
      chem3dformat
      chemdrawct
      chemtoolformat
      cifformat
      crkformat
      cssrformat
      fastsearchformat
      fastaformat
      featformat
      fhformat
      fingerprintformat
      freefracformat
      ghemicalformat
      gromos96format
      mmcifformat
      mmodformat
      moldenformat
      mpdformat
      mpqcformat
      msiformat
      msmsformat
      opendxformat
      outformat
      pcmodelformat
      pqrformat
      shelxformat
      thermoformat
      tinkerformat
      unichemformat
      viewmolformat
      xedformat
      xyzformat
      yasaraformat
      )

if(MSVC90 OR Boost_FOUND)
  set(formats_misc
    ${formats_misc}
    rxnformat
    chemdrawcdx
    chemkinformat
    rsmiformat
  )
endif(MSVC90 OR Boost_FOUND)

if(WITH_INCHI)
  add_definitions(-DINCHI_LINK_AS_DLL)
  if(NOT MSVC)
    add_subdirectory(inchi102)
    include_directories(${CMAKE_SOURCE_DIR}/include/inchi102)
    set(libs ${libs} inchi)
  else(NOT MSVC)
    set(libs ${libs} ${CMAKE_SOURCE_DIR}/windows-vc2005/libstdinchi.lib)
  endif(NOT MSVC)
  set(inchiformat_additional_sources getinchi.cpp)
  set(formats_common
    ${formats_common}
    inchiformat
  )
endif(WITH_INCHI)

if(LIBXML2_FOUND)
  if(NOT MSVC)
    include_directories(${LIBXML2_INCLUDE_DIR})
  endif(NOT MSVC)
  set(formats_xml
    cdxmlformat
    cmlformat
    pubchem
    xmlformat
  )
  if(MSVC90 OR Boost_FOUND)
    set(formats_xml
        ${formats_xml}
        cmlreactformat
    )
  endif(MSVC90 OR Boost_FOUND)
endif(LIBXML2_FOUND)

if(NOT WIN32)
  set(formats_compchem
    ${formats_compchem}
    gamessukformat
  )
  set(formats_misc
    ${formats_misc}
    xtcformat
  )
endif(NOT WIN32)

if(MINIMAL_BUILD)
  set(formats
    ${formats_common}
  )
else(MINIMAL_BUILD)
  set(formats
      ${formats_common}
      ${formats_utility}
      ${formats_compchem}
      ${formats_misc}
  )
endif(MINIMAL_BUILD)

if(NOT WIN32)
  set(libs ${libs} m)
endif(NOT WIN32)

if(BUILD_SHARED)
  if(WIN32)
    set(openbabel_srcs ${openbabel_srcs}
        dlhandler_win32
        )
  else(WIN32)
    set(openbabel_srcs ${openbabel_srcs}
        dlhandler_unix
        )
  endif(WIN32)
endif(BUILD_SHARED)

if(MSVC)
  foreach(formatgroup formats_common formats_utility formats_compchem formats_misc)
    set(formatsources "")
    foreach(format ${${formatgroup}})
      set(formatsources ${formatsources} ${format}.cpp ${${format}_additional_sources})
    endforeach(format)
    add_library(${formatgroup} ${PLUGIN_TYPE} ${formatsources}
        "${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h")
    target_link_libraries(${formatgroup} ${libs} ${BABEL_LIBRARY})
    install(TARGETS ${formatgroup}
                  RUNTIME DESTINATION ${BIN_INSTALL_DIR}
                  LIBRARY DESTINATION ${OB_PLUGIN_INSTALL_DIR}
                  ARCHIVE DESTINATION ${OB_PLUGIN_INSTALL_DIR})
    set_target_properties(${formatgroup} PROPERTIES
                  OUTPUT_NAME ${formatgroup}
                  PREFIX ""
                  SUFFIX ${MODULE_EXTENSION}
                  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endforeach(formatgroup)
else(MSVC)
foreach(format ${formats})
  add_library(${format} ${PLUGIN_TYPE} ${format}.cpp
    ${${format}_additional_sources}
    "${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h")
  target_link_libraries(${format} ${libs} ${BABEL_LIBRARY})
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

if(MSVC)
  set(xmlsources "")
  foreach(format ${formats_xml})
    set(xmlsources ${xmlsources} xml/${format}.cpp)
  endforeach(format)
  add_library(formats_xml ${PLUGIN_TYPE} ${xmlsources} xml/xml.cpp)
  target_link_libraries(formats_xml ${libs} ${BABEL_LIBRARY})
  install(TARGETS formats_xml
                  RUNTIME DESTINATION ${BIN_INSTALL_DIR}
                  LIBRARY DESTINATION ${OB_PLUGIN_INSTALL_DIR}
                  ARCHIVE DESTINATION ${OB_PLUGIN_INSTALL_DIR})
  set_target_properties(formats_xml PROPERTIES
                  OUTPUT_NAME formats_xml
                  PREFIX ""
                  SUFFIX ${MODULE_EXTENSION}
                  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
else(MSVC)
foreach(format ${formats_xml})
  add_library(${format} ${PLUGIN_TYPE} xml/${format}.cpp xml/xml.cpp)
  target_link_libraries(${format} ${libs} ${BABEL_LIBRARY})
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

