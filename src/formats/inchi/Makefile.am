
# Not sure if this is actually needed -GRH
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/data

lib_LTLIBRARIES = libinchi.la

libinchi_la_SOURCES = \
	IChICan2.c IChICano.c IChIPrt1.c IChIPrt2.c \
	IChIPrt3.c IchiMak2.c IchiMake.c IchiNorm.c \
	IchiParm.c IchiTaut.c RunICHI.c ichi_bns.c \
	ichicans.c ichiisot.c ichilnct.c ichimap1.c \
	ichimap2.c ichimap4.c ichiqueu.c ichiring.c \
	ichisort.c ichister.c inchi_dll.c inchi_dll_main.c \
	strutil.c util.c

libinchiincludedir = $(includedir)/inchi
# Only install the required header to expose to the outside world
libinchiinclude_HEADERS = inchi_api.h

EXTRA_DIST = README.html COPYING.LIB \
	aux2atom.h comdef.h extr_ct.h ichi.h ichi_bns.h \
	ichicano.h ichicant.h ichicomn.h ichicomp.h \
	ichidrp.h ichierr.h ichimain.h ichimake.h \
	ichinorm.h ichiparm.h ichiring.h ichisize.h \
	lreadmol.h mode.h strutil.h util.h inpdef.h \
	ichister.h ichitaut.h ichitime.h \
	inchi_dll.h inchi_dll_main.h
