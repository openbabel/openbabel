# Library versioning
set(SOVERSION 4)
set(LIBRARY_VERSION 4.0.0)

set(charges
        charges/gasteiger.cpp
        charges/mmff94.cpp
)

if (EIGEN2_FOUND)
# if eigen is available, compile QEq and QTPIE charge models
set(charges ${charges}
        charges/qeq.cpp
        charges/qtpie.cpp
)
endif(EIGEN2_FOUND)

set(depict_srcs
    depict/depict.cpp
    depict/svgpainter.cpp
)

set(descriptors
  descriptors/cansmidescriptor.cpp
  descriptors/cmpdfilter.cpp
  descriptors/groupcontrib.cpp
  descriptors/filters.cpp
  descriptors/inchidescriptor.cpp
  descriptors/smartsdescriptors.cpp
)

set(fingerprints
    fingerprints/finger2.cpp
    fingerprints/finger3.cpp
)

set(forcefields
    forcefields/forcefieldgaff.cpp
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

set(ops
  ops/addpolarh.cpp
  ops/gen2D.cpp
  ops/fillUC.cpp
  ops/gen3d.cpp
  ops/loader.cpp
  ops/optransform.cpp
  ops/partialcharges.cpp
  ops/sort.cpp
  ops/xout.cpp
)

set(stereo_srcs
    stereo/stereo.cpp
    stereo/tetranonplanar.cpp
    stereo/tetraplanar.cpp
    stereo/squareplanar.cpp
    stereo/cistrans.cpp
    stereo/tetrahedral.cpp
    stereo/perception.cpp
    stereo/facade.cpp
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
    chargemodel.cpp
    chiral.cpp
    data.cpp
    descriptor.cpp
    fingerprint.cpp
    forcefield.cpp
    format.cpp
    generic.cpp
    graphsym.cpp
    grid.cpp
    griddata.cpp
    kekulize.cpp
    locale.cpp
    matrix.cpp
    mcdlutil.cpp
    molchrg.cpp
    mol.cpp
    obconversion.cpp
    oberror.cpp
    obiter.cpp
    obutil.cpp
    op.cpp
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
  if(MSVC)
    set(openbabel_srcs ${openbabel_srcs}
        dlhandler_win32.cpp
        )
  else(MSVC)
    set(openbabel_srcs ${openbabel_srcs}
        dlhandler_unix.cpp
        )
  endif(MSVC)
endif(BUILD_SHARED)

set(openbabel_library_srcs
  ${openbabel_srcs}
  ${depict_srcs}
  ${math_srcs}
  ${stereo_srcs}
  ${headers}
  "${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h"
  )
if(NOT MSVC)
  set(openbabel_library_srcs
    ${openbabel_library_srcs}
    ${charges}
    ${descriptors}
    ${fingerprints}
    ${forcefields}
    ${ops}
    )
endif(NOT MSVC)

add_library(openbabel ${BUILD_TYPE}
  ${openbabel_library_srcs}
  )

if(NOT WIN32)
  set(libs ${libs} m dl)
endif(NOT WIN32)

if(ZLIB_FOUND)
  target_link_libraries(openbabel ${libs} ${ZLIB_LIBRARY})
else(ZLIB_FOUND) # omit zlib reference
  target_link_libraries(openbabel ${libs})
endif(ZLIB_FOUND)
set_target_properties(openbabel PROPERTIES
                                       VERSION ${LIBRARY_VERSION}
                                       SOVERSION ${SOVERSION}
                                       DEFINE_SYMBOL MAKE_OBDLL
                      )
if(MSVC)
  set_target_properties(openbabel PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif(MSVC)

install(TARGETS openbabel
  RUNTIME DESTINATION ${BIN_INSTALL_DIR}
  LIBRARY DESTINATION ${LIB_INSTALL_DIR}
  ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
)

if(MSVC)
  foreach(plugingroup descriptors fingerprints forcefields ops charges)
    add_library(plugin_${plugingroup} ${PLUGIN_TYPE} ${${plugingroup}}
      "${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h")
    target_link_libraries(plugin_${plugingroup} ${libs} openbabel)
    install(TARGETS plugin_${plugingroup}
      RUNTIME DESTINATION ${BIN_INSTALL_DIR}
      LIBRARY DESTINATION ${OB_PLUGIN_INSTALL_DIR}
      ARCHIVE DESTINATION ${OB_PLUGIN_INSTALL_DIR})
    set_target_properties(plugin_${plugingroup} PROPERTIES
      OUTPUT_NAME plugin_${plugingroup}
      PREFIX ""
      SUFFIX ${MODULE_EXTENSION}
      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endforeach(plugingroup)
endif(MSVC)

add_subdirectory(formats)
