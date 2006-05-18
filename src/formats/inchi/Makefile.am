
# Not sure if this is actually needed -GRH
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/data

lib_LTLIBRARIES = libinchi.la

libinchi_la_LDFLAGS = -no-undefined -version-info 0:1:0
libinchi_la_SOURCES = \
	ichi_bns.c ichican2.c ichicano.c ichicans.c ichiisot.c \
	ichilnct.c ichimak2.c ichimake.c ichimap1.c ichimap2.c \
	ichimap4.c ichinorm.c ichiparm.c ichiprt1.c ichiprt2.c \
	ichiprt3.c ichiqueu.c ichiread.c ichiring.c ichirvr1.c \
	ichirvr2.c ichirvr3.c ichirvr4.c ichirvr5.c ichirvr6.c \
	ichirvr7.c ichisort.c ichister.c ichitaut.c inchi_dll.c \
	inchi_dll_main.c runichi.c strutil.c util.c 

libinchiincludedir = $(includedir)/inchi
# Only install the required header to expose to the outside world
libinchiinclude_HEADERS = inchi_api.h

EXTRA_DIST = README.html COPYING.LIB \
	aux2atom.h comdef.h extr_ct.h ichi.h ichi_bns.h ichicano.h \
	ichicant.h ichicomn.h ichicomp.h ichidrp.h ichierr.h ichimain.h \
	ichimake.h ichinorm.h ichiparm.h ichiring.h ichirvrs.h ichisize.h \
	ichister.h ichitaut.h ichitime.h inchi_dll_main.h inchicmp.h inpdef.h \
	lreadmol.h mode.h strutil.h util.h

