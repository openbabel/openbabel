project(inchi)
set(INCHI_MAJ_VER   1)
set(INCHI_MIN_VER   07)
set(INCHI_PATCH_VER 4)

file(GLOB inchi_srcs "*.c")

# inchi code has lots of warnings, so let's ignore them
add_compile_options(-w)

# Vendored InChI 1.07.5 trips UBSAN on signed-shift sentinels (e.g.
# ichicano.c:133, "left shift of LONG_MAX by 1 cannot be represented in
# clock_t") — upstream timer code that openbabel does not author and
# cannot meaningfully fix in-tree. Disable UBSAN on this target only;
# ASAN remains active so real memory bugs in InChI are still caught.
# The flag is a no-op on builds without -fsanitize=undefined.
if(NOT MSVC)
  add_compile_options(-fno-sanitize=undefined)
endif()

if(NOT WIN32)
  set(libs ${libs} c m)
endif(NOT WIN32)

include_directories(${openbabel_SOURCE_DIR}/include/inchi)
add_library(inchi ${BUILD_TYPE} ${inchi_srcs})
target_compile_definitions(inchi PUBLIC -DTARGET_API_LIB)
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
