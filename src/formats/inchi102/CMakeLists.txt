project(inchi)
set(INCHI_MAJ_VER   0)
set(INCHI_MIN_VER   2)
set(INCHI_PATCH_VER 1)

file(GLOB inchi_srcs "*.c")

if(NOT WIN32)
  set(libs ${libs} c)
endif(NOT WIN32)

include_directories(${openbabel_SOURCE_DIR}/include/inchi102)
add_library(inchi SHARED ${inchi_srcs})
target_link_libraries(inchi ${libs})
set_target_properties(inchi PROPERTIES
                            VERSION ${INCHI_MAJ_VER}.${INCHI_MIN_VER}.${INCHI_PATCH_VER}
                            SOVERSION ${INCHI_MAJ_VER}
                            DEFINE_SYMBOL _USRDLL
)
install(TARGETS inchi
                RUNTIME DESTINATION ${BIN_INSTALL_DIR}
                LIBRARY DESTINATION ${LIB_INSTALL_DIR}
                ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
)
