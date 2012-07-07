# Library versioning
set(SOVERSION 4)
set(LIBRARY_VERSION 4.0.2)

set(charges
  charges/eem.cpp
  charges/gasteiger.cpp
  charges/none.cpp
  charges/mmff94.cpp
)

if (EIGEN2_FOUND OR EIGEN3_FOUND)
# if eigen is available, compile QEq and QTPIE charge models
set(charges ${charges}
  charges/qeq.cpp
  charges/qtpie.cpp
)
endif()

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

if(EIGEN2_FOUND OR EIGEN3_FOUND)
  set(math_srcs ${math_srcs}
    math/align.cpp
  )
endif()

set(ops
  ops/addfilename.cpp
  ops/addinindex.cpp
  ops/addpolarh.cpp
  ops/canonical.cpp
  ops/gen2D.cpp
  ops/fillUC.cpp
  ops/forcefield.cpp
  ops/gen3d.cpp
  ops/largest.cpp
  ops/loader.cpp
  ops/optransform.cpp
  ops/partialcharges.cpp
  ops/readconformers.cpp
  ops/sort.cpp
  ops/opisomorph.cpp
  ops/ophighlight.cpp
  ops/xout.cpp
)
if(EIGEN2_FOUND OR EIGEN3_FOUND)
  set(ops ${ops}
    ops/conformer.cpp
    ops/opalign.cpp
  )
endif()

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
  isomorphism.cpp
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
  query.cpp
  rand.cpp
  residue.cpp
  ring.cpp
  rotamer.cpp
  rotor.cpp
  spectrophore.cpp
  tokenst.cpp
  transform.cpp
  typer.cpp
  obmolecformat.cpp
)
if(EIGEN2_FOUND OR EIGEN3_FOUND)
  set(openbabel_srcs ${openbabel_srcs}
    conformersearch.cpp
  )
endif ()

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
else(BUILD_SHARED)
  include(formats/formats.cmake)
  foreach(format ${formats})
    set(openbabel_srcs ${openbabel_srcs} formats/${format}.cpp)
  endforeach(format ${formats})

  if(LIBXML2_FOUND AND WITH_STATIC_LIBXML)
    set(openbabel_srcs ${openbabel_srcs} formats/xml/xml.cpp)
    foreach(format ${formats_xml})
      set(openbabel_srcs ${openbabel_srcs} formats/xml/${format}.cpp)
    endforeach(format ${formats_xml})
  endif(LIBXML2_FOUND AND WITH_STATIC_LIBXML)

  foreach(plugingroup descriptors fingerprints forcefields ops charges)
    set(openbabel_srcs ${openbabel_srcs} ${${plugingroup}})
  endforeach(plugingroup descriptors fingerprints forcefields ops charges)

  if(WITH_STATIC_INCHI)
    add_subdirectory(formats/inchi103)
    set(openbabel_srcs ${openbabel_srcs} formats/getinchi.cpp ops/unique.cpp)
    set(libs ${libs} inchi)
  endif(WITH_STATIC_INCHI)
endif(BUILD_SHARED)

set(openbabel_library_srcs
  ${openbabel_srcs}
  ${depict_srcs}
  ${math_srcs}
  ${stereo_srcs}
  ${headers}
  "${openbabel_BINARY_DIR}/include/openbabel/babelconfig.h"
  )

if(WIN32)
  set(libs ${libs} ${ZLIB_LIBRARY})
else(WIN32)
  # the C math library
  if(BUILD_MIXED)
    set(libs ${libs} m)
  else(BUILD_MIXED)
    # Make sure to find the static library with .a suffix (see CMAKE_FIND_LIBRARY_SUFFIXES)
    find_library(MATH_LIBRARY NAMES m)
    if(MATH_LIBRARY)
      set(libs ${libs} ${MATH_LIBRARY})
    endif(MATH_LIBRARY)
  endif(BUILD_MIXED)

  if(BUILD_SHARED)
    # dl library if it is found and we are doing a shared build
    if(HAVE_LIBDL)
      set(libs ${libs} dl)
    endif(HAVE_LIBDL)
  endif(BUILD_SHARED)

  if(ZLIB_FOUND)
    set(libs ${libs} ${ZLIB_LIBRARY})
  endif(ZLIB_FOUND)
endif(WIN32)

add_library(openbabel ${BUILD_TYPE}
  ${openbabel_library_srcs}
  )

target_link_libraries(openbabel ${libs})

set_target_properties(openbabel PROPERTIES
  VERSION ${LIBRARY_VERSION}
  SOVERSION ${SOVERSION}
  DEFINE_SYMBOL MAKE_OBDLL
  )
if(MSVC)
  # Set the output name to openbabel-2 (preserves old behaviour). Is this desired?
  set_target_properties(openbabel PROPERTIES
    OUTPUT_NAME openbabel-2)
endif(MSVC)

# Export the openbabel target, for use in other projects
export(TARGETS openbabel APPEND FILE "${OB_EXPORTS_FILE}")

install(TARGETS openbabel
  EXPORT openbabel-targets
  RUNTIME DESTINATION ${BIN_INSTALL_DIR}
  LIBRARY DESTINATION ${LIB_INSTALL_DIR}
  ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
)

if(BUILD_SHARED)
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
      SUFFIX ${MODULE_EXTENSION})
  endforeach(plugingroup)

  add_subdirectory(formats)
endif(BUILD_SHARED)
