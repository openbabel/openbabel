// file formats
#define WITH_COMMON_FORMATS
#define WITH_UTILITY_FORMATS
//#define WITH_COMPCHEM_FORMATS
//#define WITH_MISC_FORMATS
// other plugins
#define WITH_OPS
#define WITH_DESCRIPTORS
#define WITH_FINGERPRINTS
#define WITH_FORCEFIELDS



#ifdef WITH_COMMON_FORMATS
  #include "../src/formats/smilesformat.cpp"
  #include "../src/formats/mdlformat.cpp"
  #include "../src/formats/mol2format.cpp"
  #include "../src/formats/pdbformat.cpp"
#endif

#ifdef WITH_UTILITY_FORMATS
  #include "../src/formats/copyformat.cpp"
  #include "../src/formats/molreport.cpp"
  #include "../src/formats/povrayformat.cpp"
  #include "../src/formats/reportformat.cpp"
  #include "../src/formats/titleformat.cpp"
#endif

#ifdef WITH_COMPCHEM_FORMATS
  #include "../src/formats/adfformat.cpp"
  #include "../src/formats/dmolformat.cpp"
  #include "../src/formats/fchkformat.cpp"
  #include "../src/formats/gamessformat.cpp"
  #include "../src/formats/gaussformat.cpp"
  #include "../src/formats/gausscubeformat.cpp"
  #include "../src/formats/gausszmatformat.cpp"
  #include "../src/formats/hinformat.cpp"
  #include "../src/formats/jaguarformat.cpp"
  #include "../src/formats/mopacformat.cpp"
  #include "../src/formats/nwchemformat.cpp"
  #include "../src/formats/qchemformat.cpp"
  #include "../src/formats/turbomoleformat.cpp"
  #include "../src/formats/zindoformat.cpp"
#endif

#ifdef WITH_MISC_FORMATS
  #include "../src/formats/APIInterface.cpp"
  #include "../src/formats/CSRformat.cpp"
  #include "../src/formats/PQSformat.cpp"
  #include "../src/formats/MCDLformat.cpp"
  #include "../src/formats/alchemyformat.cpp"
  #include "../src/formats/acrformat.cpp"
  #include "../src/formats/amberformat.cpp"
  #include "../src/formats/balstformat.cpp"
  #include "../src/formats/bgfformat.cpp"
  #include "../src/formats/boxformat.cpp"
  #include "../src/formats/cacaoformat.cpp"
  #include "../src/formats/cacheformat.cpp"
  #include "../src/formats/carformat.cpp"
  #include "../src/formats/cccformat.cpp"
  #include "../src/formats/chem3dformat.cpp"
  #include "../src/formats/chemdrawct.cpp"
  #include "../src/formats/chemtoolformat.cpp"
  #include "../src/formats/cifformat.cpp"
  #include "../src/formats/crkformat.cpp"
  #include "../src/formats/cssrformat.cpp"
  #include "../src/formats/fastsearchformat.cpp"
  #include "../src/formats/fastaformat.cpp"
  #include "../src/formats/featformat.cpp"
  #include "../src/formats/fhformat.cpp"
  #include "../src/formats/fingerprintformat.cpp"
  #include "../src/formats/freefracformat.cpp"
  #include "../src/formats/ghemicalformat.cpp"
  #include "../src/formats/gromos96format.cpp"
  #include "../src/formats/mmcifformat.cpp"
  #include "../src/formats/mmodformat.cpp"
  #include "../src/formats/moldenformat.cpp"
  #include "../src/formats/mpdformat.cpp"
  #include "../src/formats/mpqcformat.cpp"
  #include "../src/formats/msiformat.cpp"
  #include "../src/formats/msmsformat.cpp"
  #include "../src/formats/opendxformat.cpp"
  #include "../src/formats/outformat.cpp"
  #include "../src/formats/pcmodelformat.cpp"
  #include "../src/formats/pqrformat.cpp"
  #include "../src/formats/shelxformat.cpp"
  #include "../src/formats/thermoformat.cpp"
  #include "../src/formats/tinkerformat.cpp"
  #include "../src/formats/unichemformat.cpp"
  #include "../src/formats/viewmolformat.cpp"
  #include "../src/formats/xedformat.cpp"
  #include "../src/formats/xyzformat.cpp"
  #include "../src/formats/yasaraformat.cpp"
#endif

#ifdef WITH_OPS
  #include "../src/ops/addpolarh.cpp"
  #include "../src/ops/gen3d.cpp"
  #include "../src/ops/loader.cpp"
  #include "../src/ops/optransform.cpp"
#endif

#ifdef WITH_DESCRIPTORS
  #include "../src/descriptors/cmpdfilter.cpp"
  #include "../src/descriptors/groupcontrib.cpp"
  #include "../src/descriptors/filters.cpp"
  #include "../src/descriptors/smartsdescriptors.cpp"
#endif

#ifdef WITH_FINGERPRINTS
  #include "../src/fingerprints/finger2.cpp"
  #include "../src/fingerprints/finger3.cpp"
#endif

#ifdef WITH_FORCEFIELDS
  #include "../src/forcefields/forcefieldghemical.cpp"
  #include "../src/forcefields/forcefieldmmff94.cpp"
  #include "../src/forcefields/forcefielduff.cpp"
#endif



/*

if(MSVC OR HAVE_REGEX_H)
  set(formats_compchem
      ${formats_compchem} gamessukformat
  )
endif(MSVC OR HAVE_REGEX_H)

if(MSVC90 OR Boost_FOUND)
  set(formats_misc
    ${formats_misc}
    rxnformat
    chemdrawcdx
    chemkinformat
    rsmiformat
  )
endif(MSVC90 OR Boost_FOUND)

set(ADD_INCHI_FORMAT FALSE)
if(WITH_INCHI)
  if (NOT MSVC)
    set(ADD_INCHI_FORMAT TRUE)
  else (NOT MSVC)
    if (EXISTS ${INCHI_LIBRARY})
      set(ADD_INCHI_FORMAT TRUE)
    else (EXISTS ${INCHI_LIBRARY})
      message("WARNING: INCHI_LIBRARY not set, or does not exist.\n....InChI format will NOT be compiled.")
    endif (EXISTS ${INCHI_LIBRARY})
  endif (NOT MSVC)
endif(WITH_INCHI)

if(ADD_INCHI_FORMAT)
  add_definitions(-DINCHI_LINK_AS_DLL)
  if(NOT MSVC)
    add_subdirectory(inchi)
    include_directories(${CMAKE_SOURCE_DIR}/include/inchi)
    set(libs ${libs} inchi)
  else(NOT MSVC)
    set(libs ${libs} ${INCHI_LIBRARY})
  endif(NOT MSVC)
  set(inchiformat_additional_sources getinchi.cpp)
  set(formats_common
    ${formats_common}
    inchiformat
  )
endif(ADD_INCHI_FORMAT)

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
  if(MSVC90 OR Boost_FOUND)
    set(formats_xml
        ${formats_xml}
        cmlreactformat
    )
  endif(MSVC90 OR Boost_FOUND)
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

*/
