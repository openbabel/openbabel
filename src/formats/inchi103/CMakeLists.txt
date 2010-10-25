project(inchi)
set(INCHI_MAJ_VER   0)
set(INCHI_MIN_VER   3)
set(INCHI_PATCH_VER 1)

file(GLOB inchi_srcs "*.c")

if(NOT WIN32)
  set(libs ${libs} c)
endif(NOT WIN32)

include_directories(${openbabel_SOURCE_DIR}/include/inchi103)
add_library(inchi ${BUILD_TYPE} ${inchi_srcs})
target_link_libraries(inchi ${libs})
set_target_properties(inchi PROPERTIES
                            VERSION ${INCHI_MAJ_VER}.${INCHI_MIN_VER}.${INCHI_PATCH_VER}
                            SOVERSION ${INCHI_MAJ_VER}
                            DEFINE_SYMBOL _USRDLL
)
export(TARGETS inchi APPEND FILE "${OB_EXPORTS_FILE}")

install(TARGETS inchi
                EXPORT openbabel-targets
                RUNTIME DESTINATION ${BIN_INSTALL_DIR}
                LIBRARY DESTINATION ${LIB_INSTALL_DIR}
                ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
)
