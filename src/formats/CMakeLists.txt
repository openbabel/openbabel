
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
  povrayformat
  reportformat
  svgformat
  textformat
  titleformat
  )
set(formats_compchem
      abinitformat
      adfformat
      castepformat
      dmolformat
      fchkformat
      fhiaimsformat
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
      pwscfformat
      qchemformat
      turbomoleformat
      vaspformat
      xsfformat
      zindoformat
  )

if(MSVC OR HAVE_REGEX_H)
  set(formats_compchem
      ${formats_compchem} gamessukformat
  )
endif(MSVC OR HAVE_REGEX_H)

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
      dlpolyformat
      fastsearchformat
      fastaformat
      featformat
      fhformat
      fingerprintformat
      freefracformat
      ghemicalformat
      gromos96format
      groformat
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
      pdbqtformat
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

if(MSVC90 OR SHARED_POINTER)
  set(formats_misc
    ${formats_misc}
    rxnformat
    chemdrawcdx
    chemkinformat
    rsmiformat
  )
endif(MSVC90 OR SHARED_POINTER)

set(optional_formatgroups "")
if(CAIRO_FOUND)
  # Cairo can generate several formats (e.g. PDF, PNG): if implemented, they
  # can be added here
  set(formats_cairo
    png2format
  )
  set(png2format_additional_sources ../depict/cairopainter.cpp)
  set(optional_formatgroups
    ${optional_formatgroups} formats_cairo
    )
  include_directories(${CAIRO_INCLUDE_DIRS})
  set(libs ${libs} ${CAIRO_LIBRARIES})
endif(CAIRO_FOUND)

set(ADD_INCHI_FORMAT FALSE)
if(WITH_INCHI)
  if(NOT OPENBABEL_USE_SYSTEM_INCHI)
    set(ADD_INCHI_FORMAT TRUE)
  else()
    if (EXISTS ${INCHI_LIBRARY})
      set(ADD_INCHI_FORMAT TRUE)
    else (EXISTS ${INCHI_LIBRARY})
      message("WARNING: INCHI_LIBRARY not set, or does not exist.\n....InChI format will NOT be compiled.")
    endif()
  endif()
endif()

if(ADD_INCHI_FORMAT)
  add_definitions(-DINCHI_LINK_AS_DLL)
  if(NOT OPENBABEL_USE_SYSTEM_INCHI)
    add_subdirectory(inchi103)
    include_directories(${CMAKE_SOURCE_DIR}/include/inchi103)
    set(libs ${libs} inchi)
  else()
    include_directories(${INCHI_INCLUDE_DIR})
    set(libs ${libs} ${INCHI_LIBRARY})
  endif()
  set(inchiformat_additional_sources getinchi.cpp ../ops/unique.cpp)
  set(formats_common
    ${formats_common}
    inchiformat
  )
endif()

if(ZLIB_FOUND)
 set(formats_utility
 ${formats_utility}
  pngformat
)
endif(ZLIB_FOUND)

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
  if(MSVC90 OR SHARED_POINTER)
    set(formats_xml
        ${formats_xml}
        cmlreactformat
    )
  endif(MSVC90 OR SHARED_POINTER)
endif(LIBXML2_FOUND)

if(HAVE_RPC_XDR_H)
  set(formats_misc
    ${formats_misc}
    xtcformat
  )
endif(HAVE_RPC_XDR_H)

if(MINIMAL_BUILD)
  set(formats
    ${formats_common}
  )
else(MINIMAL_BUILD)
  set(formats "")
  foreach(formatgroup formats_common formats_utility formats_compchem formats_misc ${optional_formatgroups})
    set(formats
        ${formats} ${${formatgroup}}
    )
  endforeach(formatgroup)
endif(MINIMAL_BUILD)

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
endif(MSVC)
endif(NOT MINIMAL_BUILD)
