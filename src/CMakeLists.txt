# Library versioning
set(SOVERSION 4)
set(LIBRARY_VERSION 4.0.0)

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
  ops/sort.cpp
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
	   ${descriptors}
	   ${fingerprints}
           ${forcefields}
           ${ops}
	   )
endif(NOT MSVC)
add_library(${BABEL_LIBRARY} ${BUILD_TYPE}
                ${openbabel_library_srcs}
)

if(NOT WIN32)
  set(libs ${libs} m dl)
endif(NOT WIN32)

target_link_libraries(${BABEL_LIBRARY} ${libs} ${ZLIB_LIBRARY})
set_target_properties(${BABEL_LIBRARY} PROPERTIES
                                       VERSION ${LIBRARY_VERSION}
                                       SOVERSION ${SOVERSION}
                                       DEFINE_SYMBOL MAKE_OBDLL
                      )
if(MSVC)
    set_target_properties(${BABEL_LIBRARY} PROPERTIES
                  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif(MSVC)

install(TARGETS ${BABEL_LIBRARY}
                RUNTIME DESTINATION ${BIN_INSTALL_DIR}
                LIBRARY DESTINATION ${LIB_INSTALL_DIR}
                ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
)

if(MSVC)
  foreach(plugingroup descriptors fingerprints forcefields ops)
    add_library(plugin_${plugingroup} ${PLUGIN_TYPE} ${${plugingroup}}
        "${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h")
    target_link_libraries(plugin_${plugingroup} ${libs} ${BABEL_LIBRARY})
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
