
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/data

noinst_LTLIBRARIES = libinchi.la
libinchi_la_SOURCES = \
	ichi_bns.c ichican2.c ichicano.c ichicans.c ichiisot.c \
	ichilnct.c ichimak2.c ichimake.c ichimap1.c ichimap2.c \
	ichimap4.c ichinorm.c ichiparm.c ichiprt1.c ichiprt2.c \
	ichiprt3.c ichiqueu.c ichiring.c ichisort.c ichister.c \
	ichitaut.c inchi_dll.c inchi_dll_main.c runichi.c \
	strutil.c util.c

libinchiincludedir = $(includedir)/openbabel/inchi
libinchiinclude_HEADERS = \
	aux2atom.h comdef.h extr_ct.h ichi.h ichi_bns.h \
	ichicano.h ichicant.h ichicomn.h ichicomp.h \
	ichidrp.h ichierr.h ichimain.h ichimake.h \
	ichinorm.h ichiparm.h ichiring.h ichisize.h \
	ichister.h ichitaut.h ichitime.h inchi_api.h \
	inchi_dll.h inchi_dll_main.h inpdef.h lreadmol.h \
	mode.h strutil.h util.h 

EXTRA_DIST = InChI-License.html

