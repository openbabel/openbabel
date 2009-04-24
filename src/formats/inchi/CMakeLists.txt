project(inchi)
set(INCHI_MAJ_VER   0)
set(INCHI_MIN_VER   0)
set(INCHI_PATCH_VER 1)

set(inchi_srcs
  ichi_bns.c ichican2.c ichicano.c ichicans.c ichiisot.c
  ichilnct.c ichimak2.c ichimake.c ichimap1.c ichimap2.c
  ichimap4.c ichinorm.c ichiparm.c ichiprt1.c ichiprt2.c
  ichiprt3.c ichiqueu.c ichiread.c ichiring.c ichirvr1.c
  ichirvr2.c ichirvr3.c ichirvr4.c ichirvr5.c ichirvr6.c
  ichirvr7.c ichisort.c ichister.c ichitaut.c inchi_dll.c
  inchi_dll_main.c runichi.c strutil.c util.c
  ikey_dll.c ikey_base26.c sha2.c
)

if(NOT WIN32)
  set(libs ${libs} c)
endif(NOT WIN32)

include_directories(${openbabel_SOURCE_DIR}/include/inchi)
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
