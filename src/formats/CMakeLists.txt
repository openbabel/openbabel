if(WITH_INCHI)
  add_definitions(-DINCHI_LINK_AS_DLL)
  add_subdirectory(inchi)
  set(formats ${formats} inchiformat )
  set(libs ${libs} inchi)
  set(inchiformat_additional_sources getinchi.cpp)
endif(WITH_INCHI)

set(formats
    APIInterface
    CSRformat
    PQSformat
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
    chemdrawcdx
    chemtoolformat
    cifformat
    copyformat
    crkformat
    cssrformat
    dmolformat
    fastsearchformat
    fastaformat
    fchkformat
    featformat
    fhformat
    fingerprintformat
    freefracformat
    gamessformat
    gaussformat
    ghemicalformat
    gromos96format
    hinformat
    jaguarformat
    mdlformat
    mmodformat
    mol2format
    molreport
    mopacformat
    mpdformat
    mpqcformat
    nwchemformat
    pcmodelformat
    pdbformat
    povrayformat
    qchemformat
    reportformat
    shelxformat
    smilesformat
    thermoformat
    tinkerformat
    titleformat
    turbomoleformat
    unichemformat
    viewmolformat
    xedformat
    xyzformat
    yasaraformat
    zindoformat
)
if(BOOST_FOUND)
  # FIXME: USE_BOOST & boost check
  set(formats ${formats} rxnformat )
endif(BOOST_FOUND)

if(LIBXML2_FOUND)
  include_directories(${LIBXML2_INCLUDE_DIR})
  set(xml_formats
    cdxmlformat
    cmlformat
    pubchem
    xmlformat
  )
  if(BOOST_FOUND)
    # FIXME: USE_BOOST & boost check
    set(xml_formats
        ${xml_formats}
        cmlreactformat
    )
  endif(BOOST_FOUND)
endif(LIBXML2_FOUND)

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
  add_definitions( -DMAKE_OBDLL )
endif(BUILD_SHARED)

foreach(format ${formats})
  add_library(${format} ${BUILD_TYPE} ${format}.cpp obmolecformat.cpp)
  target_link_libraries(${format} ${libs} ${BABEL_LIBRARY})
  install(TARGETS ${format}
                  RUNTIME DESTINATION bin
                  LIBRARY DESTINATION ${LIB}/openbabel
                  ARCHIVE DESTINATION ${LIB}/openbabel)
  set_target_properties(${format} PROPERTIES OUTPUT_NAME ${format} PREFIX "")
endforeach(format)

foreach(format ${xml_formats})
  add_library(${format} ${BUILD_TYPE} xml/${format}.cpp xml/xml.cpp obmolecformat.cpp)
  target_link_libraries(${format} ${libs} ${BABEL_LIBRARY})
  install(TARGETS ${format}
                  RUNTIME DESTINATION bin
                  LIBRARY DESTINATION ${LIB}/openbabel
                  ARCHIVE DESTINATION ${LIB}/openbabel
  )
  set_target_properties(${format} PROPERTIES OUTPUT_NAME ${format} PREFIX "")
endforeach(format)
