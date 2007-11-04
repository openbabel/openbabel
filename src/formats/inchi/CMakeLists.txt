set(inchi_srcs
	ichi_bns.c ichican2.c ichicano.c ichicans.c ichiisot.c
	ichilnct.c ichimak2.c ichimake.c ichimap1.c ichimap2.c
	ichimap4.c ichinorm.c ichiparm.c ichiprt1.c ichiprt2.c
	ichiprt3.c ichiqueu.c ichiread.c ichiring.c ichirvr1.c
	ichirvr2.c ichirvr3.c ichirvr4.c ichirvr5.c ichirvr6.c
	ichirvr7.c ichisort.c ichister.c ichitaut.c inchi_dll.c
	inchi_dll_main.c runichi.c strutil.c util.c
)
include_directories(${CMAKE_SOURCE_DIR}/include/inchi)
add_library(inchi SHARED ${inchi_srcs})
set_target_properties(inchi PROPERTIES DEFINE_SYMBOL _USRDLL)
install(TARGETS inchi
                inchi DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
)
