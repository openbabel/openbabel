
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
      acesformat
      adfformat
      aoforceformat
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
      pointcloudformat
      pwscfformat
      qchemformat
      turbomoleformat
      vaspformat
      zindoformat
  )

if(MSVC OR HAVE_REGEX_H)
  set(formats_compchem
      ${formats_compchem} gamessukformat orcaformat
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
      exyzformat
      fastsearchformat
      fastaformat
      featformat
      fhformat
      fingerprintformat
      freefracformat
      ghemicalformat
      groformat
      gromos96format
      lmpdatformat
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
      stlformat
      thermoformat
      tinkerformat
      unichemformat
      viewmolformat
      xedformat
      xsfformat
      xyzformat
      yasaraformat
      genbankformat
      )

#if(MSVC90 OR Boost_FOUND)
  set(formats_misc
    ${formats_misc}
    rxnformat
    chemdrawcdx
    chemkinformat
    rsmiformat
  )
#endif(MSVC90 OR Boost_FOUND)

if(ZLIB_FOUND)
 set(formats_utility
 ${formats_utility}
  pngformat
)
endif(ZLIB_FOUND)

if(LIBXML2_FOUND AND WITH_STATIC_LIBXML)
  if(NOT MSVC)
    include_directories(${LIBXML2_INCLUDE_DIR})
  endif(NOT MSVC)
  set(formats_xml
    cdxmlformat
    cmlformat
    pubchem
    xmlformat
  )
#  if(MSVC90 OR Boost_FOUND)
    set(formats_xml
        ${formats_xml}
        cmlreactformat
    )
#  endif(MSVC90 OR Boost_FOUND)
endif(LIBXML2_FOUND AND WITH_STATIC_LIBXML)

if(HAVE_RPC_XDR_H)
  set(formats_misc
    ${formats_misc}
    xtcformat
  )
endif(HAVE_RPC_XDR_H)

if(WITH_STATIC_INCHI)
#  add_definitions(-DINCHI_LINK_AS_DLL)
  if(NOT MSVC AND NOT OPENBABEL_USE_SYSTEM_INCHI)
    include_directories(${CMAKE_SOURCE_DIR}/include/inchi)
  endif()
  set(formats_common
    ${formats_common}
    inchiformat
  )
endif(WITH_STATIC_INCHI)


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


