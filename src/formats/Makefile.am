
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/data

EXTRA_DIST = exampleformat.cpp cmlppformat.cpp

if !BUILD_SHARED

noinst_LTLIBRARIES = libformats.la
libformats_la_SOURCES = \
	APIInterface.cpp \
	CSRformat.cpp PQSformat.cpp alchemyformat.cpp \
	amberformat.cpp balstformat.cpp bgfformat.cpp boxformat.cpp \
	cacaoformat.cpp cacheformat.cpp carformat.cpp cccformat.cpp \
	chem3dformat.cpp chemdrawformat.cpp chemtoolformat.cpp \
	cml.cpp cmlformat.cpp crkformat.cpp cssrformat.cpp \
	dmolformat.cpp fastsearchformat.cpp featformat.cpp fhformat.cpp \
	freefracformat.cpp gamessformat.cpp gaussformat.cpp 
	ghemicalformat.cpp gromos96format.cpp hinformat.cpp \
	jaguarformat.cpp mdlformat.cpp mmodformat.cpp \
	mol2format.cpp mopacformat.cpp mpqcformat.cpp nwchemformat.cpp \
	pdbformat.cpp povrayformat.cpp qchemformat.cpp reportformat.cpp \
	rxnformat.cpp shelxformat.cpp smilesformat.cpp tinkerformat.cpp \
	turbomoleformat.cpp unichemformat.cpp viewmolformat.cpp \
	xedformat.cpp xyzformat.cpp zindoformat.cpp

else BUILD_SHARED

pkglib_LTLIBRARIES = \
	APIInterface.la \
	CSRformat.la PQSformat.la alchemyformat.la \
	amberformat.la balstformat.la bgfformat.la boxformat.la \
	cacaoformat.la cacheformat.la carformat.la cccformat.la \
	chem3dformat.la chemdrawformat.la chemtoolformat.la \
	cmlformat.la crkformat.la cssrformat.la \
	dmolformat.la fastsearchformat.la featformat.la fhformat.la \
	freefracformat.la gamessformat.la gaussformat.la ghemicalformat.la \
	gromos96format.la hinformat.la \
	jaguarformat.la mdlformat.la mmodformat.la \
	mol2format.la mopacformat.la mpqcformat.la nwchemformat.la \
	pdbformat.la povrayformat.la qchemformat.la reportformat.la \
	rxnformat.la shelxformat.la smilesformat.la tinkerformat.la \
	turbomoleformat.la unichemformat.la viewmolformat.la \
	xedformat.la xyzformat.la zindoformat.la

APIInterface_la_SOURCES = APIInterface.cpp
APIInterface_la_LDFLAGS = -module -avoid-version

CSRformat_la_SOURCES = CSRformat.cpp
CSRformat_la_LDFLAGS = -module -avoid-version

PQSformat_la_SOURCES = PQSformat.cpp
PQSformat_la_LDFLAGS = -module -avoid-version 

alchemyformat_la_SOURCES = alchemyformat.cpp
alchemyformat_la_LDFLAGS = -module -avoid-version 

amberformat_la_SOURCES = amberformat.cpp
amberformat_la_LDFLAGS = -module -avoid-version 

balstformat_la_SOURCES = balstformat.cpp
balstformat_la_LDFLAGS = -module -avoid-version 

bgfformat_la_SOURCES = bgfformat.cpp
bgfformat_la_LDFLAGS = -module -avoid-version 

boxformat_la_SOURCES = boxformat.cpp
boxformat_la_LDFLAGS = -module -avoid-version 

cacaoformat_la_SOURCES = cacaoformat.cpp
cacaoformat_la_LDFLAGS = -module -avoid-version 

cacheformat_la_SOURCES = cacheformat.cpp
cacheformat_la_LDFLAGS = -module -avoid-version 

carformat_la_SOURCES = carformat.cpp
carformat_la_LDFLAGS = -module -avoid-version 

cccformat_la_SOURCES = cccformat.cpp
cccformat_la_LDFLAGS = -module -avoid-version 

chem3dformat_la_SOURCES = chem3dformat.cpp
chem3dformat_la_LDFLAGS = -module -avoid-version 

chemdrawformat_la_SOURCES = chemdrawformat.cpp
chemdrawformat_la_LDFLAGS = -module -avoid-version 

chemtoolformat_la_SOURCES = chemtoolformat.cpp
chemtoolformat_la_LDFLAGS = -module -avoid-version 

cmlformat_la_SOURCES = cml.cpp cmlformat.cpp
cmlformat_la_LDFLAGS = -module -avoid-version 

crkformat_la_SOURCES = crkformat.cpp
crkformat_la_LDFLAGS = -module -avoid-version 

cssrformat_la_SOURCES = cssrformat.cpp
cssrformat_la_LDFLAGS = -module -avoid-version 

dmolformat_la_SOURCES = dmolformat.cpp
dmolformat_la_LDFLAGS = -module -avoid-version 

fastsearchformat_la_SOURCES = fastsearchformat.cpp
fastsearchformat_la_LDFLAGS = -module -avoid-version 

featformat_la_SOURCES = featformat.cpp
featformat_la_LDFLAGS = -module -avoid-version 

fhformat_la_SOURCES = fhformat.cpp
fhformat_la_LDFLAGS = -module -avoid-version

freefracformat_la_SOURCES = freefracformat.cpp
freefracformat_la_LDFLAGS = -module -avoid-version 

gamessformat_la_SOURCES = gamessformat.cpp
gamessformat_la_LDFLAGS = -module -avoid-version 

gaussformat_la_SOURCES = gaussformat.cpp
gaussformat_la_LDFLAGS = -module -avoid-version 

ghemicalformat_la_SOURCES = ghemicalformat.cpp
ghemicalformat_la_LDFLAGS = -module -avoid-version 

gromos96format_la_SOURCES = gromos96format.cpp
gromos96format_la_LDFLAGS = -module -avoid-version 

hinformat_la_SOURCES = hinformat.cpp
hinformat_la_LDFLAGS = -module -avoid-version 

jaguarformat_la_SOURCES = jaguarformat.cpp
jaguarformat_la_LDFLAGS = -module -avoid-version 

mdlformat_la_SOURCES = mdlformat.cpp
mdlformat_la_LDFLAGS = -module -avoid-version 

#mm3format_la_SOURCES = mm3format.cpp
#mm3format_la_LDFLAGS = -module -avoid-version 

mmodformat_la_SOURCES = mmodformat.cpp
mmodformat_la_LDFLAGS = -module -avoid-version 

mol2format_la_SOURCES = mol2format.cpp
mol2format_la_LDFLAGS = -module -avoid-version 

mopacformat_la_SOURCES = mopacformat.cpp
mopacformat_la_LDFLAGS = -module -avoid-version 

mpqcformat_la_SOURCES = mpqcformat.cpp
mpqcformat_la_LDFLAGS = -module -avoid-version 

nwchemformat_la_SOURCES = nwchemformat.cpp
nwchemformat_la_LDFLAGS = -module -avoid-version 

pdbformat_la_SOURCES = pdbformat.cpp
pdbformat_la_LDFLAGS = -module -avoid-version 

povrayformat_la_SOURCES = povrayformat.cpp
povrayformat_la_LDFLAGS = -module -avoid-version 

qchemformat_la_SOURCES = qchemformat.cpp
qchemformat_la_LDFLAGS = -module -avoid-version 

reportformat_la_SOURCES = reportformat.cpp
reportformat_la_LDFLAGS = -module -avoid-version 

rxnformat_la_SOURCES = rxnformat.cpp
rxnformat_la_LDFLAGS = -module -avoid-version 

shelxformat_la_SOURCES = shelxformat.cpp
shelxformat_la_LDFLAGS = -module -avoid-version 

smilesformat_la_SOURCES = smilesformat.cpp
smilesformat_la_LDFLAGS = -module -avoid-version 

tinkerformat_la_SOURCES = tinkerformat.cpp
tinkerformat_la_LDFLAGS = -module -avoid-version 

turbomoleformat_la_SOURCES = turbomoleformat.cpp
turbomoleformat_la_LDFLAGS = -module -avoid-version 

unichemformat_la_SOURCES = unichemformat.cpp
unichemformat_la_LDFLAGS = -module -avoid-version 

viewmolformat_la_SOURCES = viewmolformat.cpp
viewmolformat_la_LDFLAGS = -module -avoid-version 

xedformat_la_SOURCES = xedformat.cpp
xedformat_la_LDFLAGS = -module -avoid-version 

xyzformat_la_SOURCES = xyzformat.cpp
xyzformat_la_LDFLAGS = -module -avoid-version 

zindoformat_la_SOURCES = zindoformat.cpp
zindoformat_la_LDFLAGS = -module -avoid-version 

endif

# InChI code cannot be distributed with Open Babel due to licensing issues
# (current InChI code is not released under any sort of free software or 
#  open source license. IUPAC is supposedly reconsidering the InChI license.)
#
# If you "happen" to have the inchi source code and drop it into the following
# subdirectory, things might work out OK. It might also work if you happen
# to have a pre-compiled libinchi that's discovered by the configure script.
#
# SUBDIRS = inchi

if BUILD_INCHI
libopenbabelincludedir = $(includedir)/openbabel-2.0/openbabel
libopenbabelinclude_HEADERS = inchi_api.h

if !BUILD_SHARED
libformats_la_SOURCES +=  \
		inchiformat.cpp
else

pkglib_LTLIBRARIES += \
		inchiformat.la

inchiformat_la_SOURCES = inchiformat.cpp
inchiformat_la_LDFLAGS = -module -avoid-version 

endif #build shared

endif #build inchi