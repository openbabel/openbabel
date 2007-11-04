if(WITH_INCHI)
  add_definitions(-DINCHI_LINK_AS_DLL)
  add_subdirectory(formats/inchi)
endif(WITH_INCHI)

set(fingerprints_srcs
    fingerprints/finger2.cpp
    fingerprints/finger3.cpp
)

set(forcefields_srcs
    forcefields/forcefieldghemical.cpp
)

set(formats_srcs
    formats/APIInterface.cpp
    formats/obmolecformat.cpp
    formats/CSRformat.cpp
    formats/PQSformat.cpp
    formats/alchemyformat.cpp
    formats/acrformat.cpp
    formats/amberformat.cpp
    formats/balstformat.cpp
    formats/bgfformat.cpp
    formats/boxformat.cpp
    formats/cacaoformat.cpp
    formats/cacheformat.cpp
    formats/carformat.cpp
    formats/cccformat.cpp
    formats/chem3dformat.cpp
    formats/chemdrawct.cpp
    formats/chemdrawcdx.cpp
    formats/chemtoolformat.cpp
    formats/cifformat.cpp
    formats/copyformat.cpp
    formats/crkformat.cpp
    formats/cssrformat.cpp
    formats/dmolformat.cpp
    formats/fastsearchformat.cpp
    formats/fastaformat.cpp
    formats/fchkformat.cpp
    formats/featformat.cpp
    formats/fhformat.cpp
    formats/fingerprintformat.cpp
    formats/freefracformat.cpp
    formats/gamessformat.cpp
    formats/gaussformat.cpp
    formats/ghemicalformat.cpp
    formats/gromos96format.cpp
    formats/hinformat.cpp
    formats/jaguarformat.cpp
    formats/mdlformat.cpp
    formats/mmodformat.cpp
    formats/mol2format.cpp
    formats/molreport.cpp
    formats/mopacformat.cpp
    formats/mpdformat.cpp
    formats/mpqcformat.cpp
    formats/nwchemformat.cpp
    formats/pcmodelformat.cpp
    formats/pdbformat.cpp
    formats/povrayformat.cpp
    formats/qchemformat.cpp
    formats/reportformat.cpp
    formats/rxnformat.cpp
    formats/shelxformat.cpp
    formats/smilesformat.cpp
    formats/thermoformat.cpp
    formats/tinkerformat.cpp
    formats/titleformat.cpp
    formats/turbomoleformat.cpp
    formats/unichemformat.cpp
    formats/viewmolformat.cpp
    formats/xedformat.cpp
    formats/xyzformat.cpp
    formats/yasaraformat.cpp
    formats/zindoformat.cpp
)

if(WITH_INCHI)
  set(formats_srcs
    ${formats_srcs}
	formats/getinchi.cpp
	formats/inchiformat.cpp
  )
  set(libs ${libs} inchi)
endif(WITH_INCHI)

if(LIBXML2_FOUND)
  include_directories(${LIBXML2_INCLUDE_DIR})
  set(xml_srcs
    formats/xml/cdxmlformat.cpp
    formats/xml/cmlformat.cpp
    formats/xml/cmlreactlformat.cpp
    formats/xml/pubchem.cpp
    formats/xml/xmlformat.cpp
    formats/xml/xml.cpp
  )
endif(LIBXML2_FOUND)

set(math_srcs
    math/matrix3x3.cpp
    math/vector3.cpp
)

set(openbabel_srcs
    atom.cpp
    base.cpp
    bitvec.cpp
    bond.cpp
    bondtyper.cpp
    canon.cpp
    chains.cpp
    chiral.cpp
    data.cpp
    fingerprint.cpp
    forcefield.cpp
    generic.cpp
    grid.cpp
    kekulize.cpp
    matrix.cpp
    mol.cpp
    molchrg.cpp
    obconversion.cpp
    oberror.cpp
    obiter.cpp
    obutil.cpp
    parsmart.cpp
    patty.cpp
    phmodel.cpp
    rand.cpp
    residue.cpp
    ring.cpp
    rotor.cpp
    rotamer.cpp
    tokenst.cpp
    typer.cpp
    transform.cpp
    groupcontrib.cpp
)

if(BUILD_SHARED)
  if(WIN32)
    set(openbabel_srcs ${openbabel_srcs}
        dlhandler_win32.cpp
        )
  else(WIN32)
    set(openbabel_srcs ${openbabel_srcs}
        dlhandler_unix.cpp
        )
  endif(WIN32)
endif(BUILD_SHARED)

add_library(openbabel-2 ${BUILD_TYPE}
                ${openbabel_srcs}
                ${fingerprints_srcs}
                ${forcefields_srcs}
                ${formats_srcs}
                ${xml_srcs}
                ${math_srcs}
)

if(NOT WIN32)
  set(libs ${libs} m)
endif(NOT WIN32)

target_link_libraries(openbabel-2 ${libs} ${LIBXML2_LIBRARIES} ${ZLIB_LIBRARY})
install(TARGETS openbabel-2
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
)
