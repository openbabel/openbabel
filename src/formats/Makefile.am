INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/data
# Make sure to install the inchi_api.h header
libopenbabelincludedir = $(includedir)/openbabel-2.0/openbabel
libopenbabelinclude_HEADERS = inchi_api.h obmolecformat.h

EXTRA_DIST = exampleformat.cpp
SUBDIRS =

if BUILD_XML
SUBDIRS += xml
endif

if BUILD_INCHI
SUBDIRS += inchi
endif

# Unfortunately there are a few formats which still cannot compile shared
# modules (e.g., Cygwin). So we need to keep a duplicate target for "libformats"

if !BUILD_SHARED
noinst_LTLIBRARIES = libformats.la
if BUILD_INCHI
libformats_la_LIBADD = inchi/libinchi.la
else
libformats_la_LIBADD = -linchi
endif
libformats_la_SOURCES = \
	APIInterface.cpp obmolecformat.cpp \
	CSRformat.cpp PQSformat.cpp alchemyformat.cpp \
	amberformat.cpp balstformat.cpp bgfformat.cpp boxformat.cpp \
	cacaoformat.cpp cacheformat.cpp carformat.cpp cccformat.cpp \
	chem3dformat.cpp chemdrawformat.cpp chemtoolformat.cpp \
	copyformat.cpp crkformat.cpp cssrformat.cpp \
	dmolformat.cpp fastsearchformat.cpp fastaformat.cpp \
	featformat.cpp fhformat.cpp fingerprintformat.cpp \
	freefracformat.cpp gamessformat.cpp gaussformat.cpp \
	ghemicalformat.cpp gromos96format.cpp hinformat.cpp \
	inchiformat.cpp jaguarformat.cpp mdlformat.cpp mmodformat.cpp \
	mol2format.cpp mopacformat.cpp mpdformat.cpp mpqcformat.cpp \
	nwchemformat.cpp pcmodelformat.cpp \
	pdbformat.cpp povrayformat.cpp qchemformat.cpp reportformat.cpp \
	rxnformat.cpp shelxformat.cpp smilesformat.cpp thermoformat.cpp \
	tinkerformat.cpp titleformat.cpp \
	turbomoleformat.cpp unichemformat.cpp viewmolformat.cpp \
	xedformat.cpp xyzformat.cpp yasaraformat.cpp zindoformat.cpp

else BUILD_SHARED

pkglib_LTLIBRARIES = \
	APIInterface.la \
	CSRformat.la PQSformat.la alchemyformat.la \
	amberformat.la balstformat.la bgfformat.la boxformat.la \
	cacaoformat.la cacheformat.la carformat.la cccformat.la \
	chem3dformat.la chemdrawformat.la chemtoolformat.la \
	copyformat.la crkformat.la cssrformat.la \
	dmolformat.la fastsearchformat.la fastaformat.la \
	featformat.la fhformat.la fingerprintformat.la \
	freefracformat.la gamessformat.la gaussformat.la ghemicalformat.la \
	gromos96format.la hinformat.la inchiformat.la \
	jaguarformat.la mdlformat.la mmodformat.la mpdformat.la \
	mol2format.la mopacformat.la mpqcformat.la nwchemformat.la \
	pcmodelformat.la \
	pdbformat.la povrayformat.la qchemformat.la reportformat.la \
	rxnformat.la shelxformat.la smilesformat.la thermoformat.la \
	tinkerformat.la titleformat.la \
	turbomoleformat.la unichemformat.la viewmolformat.la \
	xedformat.la xyzformat.la yasaraformat.la zindoformat.la

APIInterface_la_SOURCES = APIInterface.cpp
APIInterface_la_LDFLAGS = -module -avoid-version

CSRformat_la_SOURCES = CSRformat.cpp
CSRformat_la_LDFLAGS = -module -avoid-version

PQSformat_la_SOURCES = PQSformat.cpp obmolecformat.cpp
PQSformat_la_LDFLAGS = -module -avoid-version 

alchemyformat_la_SOURCES = alchemyformat.cpp obmolecformat.cpp
alchemyformat_la_LDFLAGS = -module -avoid-version 

amberformat_la_SOURCES = amberformat.cpp obmolecformat.cpp
amberformat_la_LDFLAGS = -module -avoid-version 

balstformat_la_SOURCES = balstformat.cpp obmolecformat.cpp
balstformat_la_LDFLAGS = -module -avoid-version 

bgfformat_la_SOURCES = bgfformat.cpp obmolecformat.cpp
bgfformat_la_LDFLAGS = -module -avoid-version 

boxformat_la_SOURCES = boxformat.cpp obmolecformat.cpp
boxformat_la_LDFLAGS = -module -avoid-version 

cacaoformat_la_SOURCES = cacaoformat.cpp obmolecformat.cpp
cacaoformat_la_LDFLAGS = -module -avoid-version 

cacheformat_la_SOURCES = cacheformat.cpp
cacheformat_la_LDFLAGS = -module -avoid-version 

carformat_la_SOURCES = carformat.cpp obmolecformat.cpp
carformat_la_LDFLAGS = -module -avoid-version 

cccformat_la_SOURCES = cccformat.cpp obmolecformat.cpp
cccformat_la_LDFLAGS = -module -avoid-version 

chem3dformat_la_SOURCES = chem3dformat.cpp obmolecformat.cpp
chem3dformat_la_LDFLAGS = -module -avoid-version 

chemdrawformat_la_SOURCES = chemdrawformat.cpp obmolecformat.cpp
chemdrawformat_la_LDFLAGS = -module -avoid-version 

chemtoolformat_la_SOURCES = chemtoolformat.cpp
chemtoolformat_la_LDFLAGS = -module -avoid-version 

copyformat_la_SOURCES = copyformat.cpp
copyformat_la_LDFLAGS = -module -avoid-version 

crkformat_la_SOURCES = crkformat.cpp obmolecformat.cpp
crkformat_la_LDFLAGS = -module -avoid-version 

cssrformat_la_SOURCES = cssrformat.cpp
cssrformat_la_LDFLAGS = -module -avoid-version 

dmolformat_la_SOURCES = dmolformat.cpp obmolecformat.cpp
dmolformat_la_LDFLAGS = -module -avoid-version 

fastsearchformat_la_SOURCES = fastsearchformat.cpp
fastsearchformat_la_LDFLAGS = -module -avoid-version 

fastaformat_la_SOURCES = fastaformat.cpp obmolecformat.cpp
fastaformat_la_LDFLAGS = -module -avoid-version 

featformat_la_SOURCES = featformat.cpp obmolecformat.cpp
featformat_la_LDFLAGS = -module -avoid-version 

fhformat_la_SOURCES = fhformat.cpp
fhformat_la_LDFLAGS = -module -avoid-version

fingerprintformat_la_SOURCES = fingerprintformat.cpp obmolecformat.cpp
fingerprintformat_la_LDFLAGS = -module -avoid-version

freefracformat_la_SOURCES = freefracformat.cpp obmolecformat.cpp
freefracformat_la_LDFLAGS = -module -avoid-version 

gamessformat_la_SOURCES = gamessformat.cpp obmolecformat.cpp
gamessformat_la_LDFLAGS = -module -avoid-version 

gaussformat_la_SOURCES = gaussformat.cpp obmolecformat.cpp
gaussformat_la_LDFLAGS = -module -avoid-version 

ghemicalformat_la_SOURCES = ghemicalformat.cpp obmolecformat.cpp
ghemicalformat_la_LDFLAGS = -module -avoid-version 

gromos96format_la_SOURCES = gromos96format.cpp
gromos96format_la_LDFLAGS = -module -avoid-version 

hinformat_la_SOURCES = hinformat.cpp obmolecformat.cpp
hinformat_la_LDFLAGS = -module -avoid-version 

inchiformat_la_SOURCES = inchiformat.cpp obmolecformat.cpp
if BUILD_INCHI
inchiformat_la_LDFLAGS = -module -avoid-version inchi/libinchi.la
else
inchiformat_la_LDFLAGS = -module -avoid-version
inchiformat_la_LIBADD = -linchi
endif

jaguarformat_la_SOURCES = jaguarformat.cpp obmolecformat.cpp
jaguarformat_la_LDFLAGS = -module -avoid-version 

mdlformat_la_SOURCES = mdlformat.cpp obmolecformat.cpp
mdlformat_la_LDFLAGS = -module -avoid-version 

mmodformat_la_SOURCES = mmodformat.cpp obmolecformat.cpp
mmodformat_la_LDFLAGS = -module -avoid-version 

mol2format_la_SOURCES = mol2format.cpp obmolecformat.cpp
mol2format_la_LDFLAGS = -module -avoid-version 

mopacformat_la_SOURCES = mopacformat.cpp obmolecformat.cpp
mopacformat_la_LDFLAGS = -module -avoid-version 

mpdformat_la_SOURCES = mpdformat.cpp obmolecformat.cpp
mpdformat_la_LDFLAGS = -module -avoid-version 

mpqcformat_la_SOURCES = mpqcformat.cpp obmolecformat.cpp
mpqcformat_la_LDFLAGS = -module -avoid-version 

nwchemformat_la_SOURCES = nwchemformat.cpp obmolecformat.cpp
nwchemformat_la_LDFLAGS = -module -avoid-version 

pcmodelformat_la_SOURCES = pcmodelformat.cpp obmolecformat.cpp
pcmodelformat_la_LDFLAGS = -module -avoid-version 

pdbformat_la_SOURCES = pdbformat.cpp obmolecformat.cpp
pdbformat_la_LDFLAGS = -module -avoid-version 

povrayformat_la_SOURCES = povrayformat.cpp
povrayformat_la_LDFLAGS = -module -avoid-version 

qchemformat_la_SOURCES = qchemformat.cpp obmolecformat.cpp
qchemformat_la_LDFLAGS = -module -avoid-version 

reportformat_la_SOURCES = reportformat.cpp
reportformat_la_LDFLAGS = -module -avoid-version 

rxnformat_la_SOURCES = rxnformat.cpp
rxnformat_la_LDFLAGS = -module -avoid-version 

shelxformat_la_SOURCES = shelxformat.cpp obmolecformat.cpp
shelxformat_la_LDFLAGS = -module -avoid-version 

smilesformat_la_SOURCES = smilesformat.cpp obmolecformat.cpp
smilesformat_la_LDFLAGS = -module -avoid-version 

thermoformat_la_SOURCES = thermoformat.cpp obmolecformat.cpp
thermoformat_la_LDFLAGS = -module -avoid-version 

tinkerformat_la_SOURCES = tinkerformat.cpp
tinkerformat_la_LDFLAGS = -module -avoid-version 

titleformat_la_SOURCES = titleformat.cpp obmolecformat.cpp
titleformat_la_LDFLAGS = -module -avoid-version 

turbomoleformat_la_SOURCES = turbomoleformat.cpp obmolecformat.cpp
turbomoleformat_la_LDFLAGS = -module -avoid-version

unichemformat_la_SOURCES = unichemformat.cpp obmolecformat.cpp
unichemformat_la_LDFLAGS = -module -avoid-version 

viewmolformat_la_SOURCES = viewmolformat.cpp obmolecformat.cpp
viewmolformat_la_LDFLAGS = -module -avoid-version 

xedformat_la_SOURCES = xedformat.cpp
xedformat_la_LDFLAGS = -module -avoid-version 

xyzformat_la_SOURCES = xyzformat.cpp obmolecformat.cpp
xyzformat_la_LDFLAGS = -module -avoid-version 

yasaraformat_la_SOURCES = yasaraformat.cpp obmolecformat.cpp
yasaraformat_la_LDFLAGS = -module -avoid-version 

zindoformat_la_SOURCES = zindoformat.cpp
zindoformat_la_LDFLAGS = -module -avoid-version 

endif
