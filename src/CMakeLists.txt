# Library versioning
set(SOVERSION 3)
set(LIBRARY_VERSION 3.0.3)

set(descriptors_srcs
  descriptors/cmpdfilter.cpp
  descriptors/groupcontrib.cpp
  descriptors/filters.cpp
  descriptors/smartsdescriptors.cpp
)

set(fingerprints_srcs
    fingerprints/finger2.cpp
    fingerprints/finger3.cpp
)

set(forcefields_srcs
    forcefields/forcefieldghemical.cpp
    forcefields/forcefieldmmff94.cpp
    forcefields/forcefielduff.cpp
)

set(math_srcs
    math/matrix3x3.cpp
    math/spacegroup.cpp
    math/transform3d.cpp
    math/vector3.cpp
)

set(ops_srcs
  ops/addpolarh.cpp
  ops/fillUC.cpp
  ops/gen3d.cpp
  ops/loader.cpp
  ops/optransform.cpp
)

set(openbabel_srcs
    alias.cpp
    atom.cpp
    base.cpp
    bitvec.cpp
    bond.cpp
    bondtyper.cpp
    builder.cpp
    canon.cpp
    chains.cpp
    chiral.cpp
    data.cpp
    descriptor.cpp
    fingerprint.cpp
    forcefield.cpp
    format.cpp
    generic.cpp
    grid.cpp
    griddata.cpp
    kekulize.cpp
    locale.cpp
    matrix.cpp
    molchrg.cpp
    mol.cpp
    obconversion.cpp
    oberror.cpp
    obiter.cpp
    obutil.cpp
    op.cpp
    ops.cpp
    parsmart.cpp
    patty.cpp
    phmodel.cpp
    plugin.cpp
    pointgroup.cpp
    rand.cpp
    residue.cpp
    ring.cpp
    rotamer.cpp
    rotor.cpp
    tokenst.cpp
    transform.cpp
    typer.cpp
    obmolecformat.cpp
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

add_library(${BABEL_LIBRARY} ${BUILD_TYPE}
                ${openbabel_srcs}
                ${descriptors_srcs}
                ${fingerprints_srcs}
                ${forcefields_srcs}
                ${math_srcs}
                ${ops_srcs}
                ${headers}
                "${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h"
)

if(NOT WIN32)
  set(libs ${libs} m dl)
endif(NOT WIN32)

target_link_libraries(${BABEL_LIBRARY} ${libs} ${LIBXML2_LIBRARIES} ${ZLIB_LIBRARY})
set_target_properties(${BABEL_LIBRARY} PROPERTIES
                                       VERSION ${LIBRARY_VERSION}
                                       SOVERSION ${SOVERSION}
                                       DEFINE_SYMBOL MAKE_OBDLL
                      )
install(TARGETS ${BABEL_LIBRARY}
                RUNTIME DESTINATION ${BIN_INSTALL_DIR}
                LIBRARY DESTINATION ${LIB_INSTALL_DIR}
                ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
)

add_subdirectory(formats)
