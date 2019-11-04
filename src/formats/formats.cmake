
set(formats_common
  smilesformat
  mdlformat
  mol2format
  pdbformat
 )

set(formats_utility
  asciiformat
  copyformat
  MNAformat
  molreport
  nulformat
  painterformat
  povrayformat
  reportformat
  svgformat
  textformat
  titleformat
  )
set(painterformat_additional_sources ../depict/commandpainter.cpp)
set(asciiformat_additional_sources   ../depict/asciipainter.cpp)
if(EIGEN2_FOUND OR EIGEN3_FOUND)
  set(formats_utility ${formats_utility}
      confabreport
     )
endif()

set(formats_compchem
      acesformat
      abinitformat
      adfformat
      aoforceformat
      castepformat
      crystal09format
      daltonformat
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
      siestaformat
      turbomoleformat
      vaspformat
      xsfformat
      zindoformat
  )

if(WITH_MAEPARSER)
    set(formats_compchem ${formats_compchem}
        maeformat
       )
endif()


if(MSVC OR HAVE_REGEX_H)
  set(formats_compchem
      ${formats_compchem} gamessukformat
  )
  set(formats_compchem
      ${formats_compchem} orcaformat
  )
endif(MSVC OR HAVE_REGEX_H)


if(WITH_JSON)
  set(formats_json
    chemdoodlejsonformat
    pubchemjsonformat
  )
endif()

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
      cofformat
      crkformat
      cssrformat
      dlpolyformat
      exyzformat      
      fastsearchformat
      fastaformat
      featformat
      fhformat
      fingerprintformat
      fpsformat
      freefracformat
      ghemicalformat
      gromos96format
      groformat
      lmpdatformat
      lpmdformat
      mdffformat
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
      pointcloudformat
      posformat
      pqrformat
      shelxformat
      smileyformat
      stlformat
      thermoformat
      tinkerformat
      unichemformat
      viewmolformat
      xedformat
      xyzformat
      yasaraformat
      )

# genbankformat can currently only be built statically
if(NOT BUILD_SHARED)
  set(formats_misc ${formats_misc} genbankformat)
endif(NOT BUILD_SHARED)

if(MSVC OR SHARED_POINTER)
  set(formats_misc
    ${formats_misc}
    rxnformat
    chemdrawcdx
    chemkinformat
    rinchiformat
    rsmiformat
  )
endif(MSVC OR SHARED_POINTER)

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

# Inchi settings for shared builds
if(BUILD_SHARED)
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
    add_definitions(-DBUILD_LINK_AS_DLL)
    if(NOT OPENBABEL_USE_SYSTEM_INCHI)
      add_subdirectory(libinchi)
      include_directories(${CMAKE_SOURCE_DIR}/include/inchi)
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

# Inchi settings for static builds
elseif(WITH_STATIC_INCHI)
  #  add_definitions(-DINCHI_LINK_AS_DLL)
  if(NOT MSVC AND NOT OPENBABEL_USE_SYSTEM_INCHI)
    include_directories(${CMAKE_SOURCE_DIR}/include/inchi)
  endif()
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

if(LIBXML2_FOUND AND (BUILD_SHARED OR WITH_STATIC_LIBXML))
  if(NOT MSVC)
    include_directories(${LIBXML2_INCLUDE_DIR})
  endif(NOT MSVC)
  set(formats_xml
    cdxmlformat
    cmlformat
    pubchem
    xmlformat
  )
  if(MSVC OR SHARED_POINTER)
    set(formats_xml
        ${formats_xml}
        cmlreactformat
    )
  endif(MSVC OR SHARED_POINTER)
endif(LIBXML2_FOUND AND (BUILD_SHARED OR WITH_STATIC_LIBXML))

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
