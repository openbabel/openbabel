INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/data -I$(top_srcdir)/include

# KFc.h and ArrayList.h for ADF formats
EXTRA_DIST = exampleformat.cpp ArrayList.h CMakeLists.txt
SUBDIRS =

if BUILD_XML
SUBDIRS += xml
endif

if BUILD_INCHI
SUBDIRS += inchi102
endif

if USE_BOOST_FLAGS
AM_CXXFLAGS = $(BOOST_CPPFLAGS)
endif

if BUILD_PCH
BUILT_SOURCES = all.h.gch
PCHFLAGS=-Winvalid-pch -x c++-header $(CPPFLAGS) $(INCLUDES)
all.h.gch: all.h Makefile \
	$(srcdir)/../../include/openbabel/babelconfig.h \
    $(srcdir)/../../include/openbabel/base.h \
	$(srcdir)/../../include/openbabel/generic.h \
	$(srcdir)/../../include/openbabel/mol.h \
	$(srcdir)/../../include/openbabel/oberror.h \
	$(srcdir)/../../include/openbabel/obconversion.h \
	$(srcdir)/../../include/openbabel/obmolecformat.h
	rm -f $@
	$(CXX) $(PCHFLAGS) $<
AM_CXXFLAGS= -include $(top_srcdir)/src/formats/all.h
clean-local:
	-rm -f all.h.gch
endif

# Unfortunately there are a few formats which still cannot compile shared
# modules (e.g., Cygwin). So we need to keep a duplicate target for "libformats"

if !BUILD_SHARED
noinst_LTLIBRARIES = libformats.la
if BUILD_INCHI
libformats_la_LIBADD = inchi102/libinchi.la
else
libformats_la_LIBADD = -linchi
endif
libformats_la_SOURCES = \
	APIInterface.cpp MCDLformat.cpp \
	CSRformat.cpp PQSformat.cpp alchemyformat.cpp \
	acrformat.cpp \
	amberformat.cpp balstformat.cpp bgfformat.cpp boxformat.cpp \
	cacaoformat.cpp cacheformat.cpp carformat.cpp cccformat.cpp \
	chem3dformat.cpp chemdrawct.cpp \
	chemkinformat.cpp \
	chemtoolformat.cpp cifformat.cpp \
	copyformat.cpp crkformat.cpp cssrformat.cpp \
	dmolformat.cpp fastsearchformat.cpp fastaformat.cpp \
	fchkformat.cpp \
	featformat.cpp fhformat.cpp fingerprintformat.cpp \
	freefracformat.cpp gamessformat.cpp gaussformat.cpp \
	gausscubeformat.cpp gausszmatformat.cpp \
	ghemicalformat.cpp gromos96format.cpp hinformat.cpp \
	getinchi.cpp \
	inchiformat.cpp jaguarformat.cpp mdlformat.cpp mmcifformat.cpp \
	mmodformat.cpp mol2format.cpp molreport.cpp mopacformat.cpp \
	moldenformat.cpp molproformat.cpp msmsformat.cpp \
	mpdformat.cpp mpqcformat.cpp nwchemformat.cpp outformat.cpp \
	pcmodelformat.cpp \
	pdbformat.cpp povrayformat.cpp pqrformat.cpp qchemformat.cpp reportformat.cpp \
	shelxformat.cpp smilesformat.cpp thermoformat.cpp \
	tinkerformat.cpp titleformat.cpp \
	turbomoleformat.cpp unichemformat.cpp viewmolformat.cpp \
	xedformat.cpp xyzformat.cpp yasaraformat.cpp zindoformat.cpp \
	rxnformat.cpp rsmiformat.cpp

else BUILD_SHARED

# First declare where these will be installed
if BUILD_SHARED_UNIX
formatsdir = $(pkglibdir)/$(VERSION)
win32ldflags = 

else
#MinGW: install with babel.exe as *.obf files
formatsdir = $(bindir)
win32ldflags = -export-symbols-regex ".*"
install-data-hook:
	(cd $(DESTDIR)$(bindir); \
	for format in APIInterface.dll *format.dll; do \
		mv $${format} `basename $${format} .dll`.obf; \
	done)
endif

formats_LTLIBRARIES = \
	APIInterface.la \
	CSRformat.la PQSformat.la MCDLformat.la alchemyformat.la \
	acrformat.la adfformat.la \
	amberformat.la balstformat.la bgfformat.la boxformat.la \
	cacaoformat.la cacheformat.la \
	carformat.la cccformat.la \
	chem3dformat.la chemdrawctformat.la \
	chemtoolformat.la cifformat.la \
	copyformat.la crkformat.la cssrformat.la \
	dmolformat.la fastsearchformat.la fastaformat.la \
	fchkformat.la \
	featformat.la fhformat.la fingerprintformat.la \
	freefracformat.la gamessformat.la \
	gaussformat.la ghemicalformat.la \
	gausscubeformat.la gausszmatformat.la \
	gromos96format.la hinformat.la \
	jaguarformat.la mdlformat.la mmcifformat.la mmodformat.la \
	mpdformat.la mol2format.la moldenformat.la molreportformat.la \
	mopacformat.la molproformat.la msiformat.la msmsformat.la \
	mpqcformat.la nwchemformat.la outformat.la \
	pcmodelformat.la pdbformat.la \
	povrayformat.la pqrformat.la opendxformat.la \
	qchemformat.la reportformat.la \
	shelxformat.la smilesformat.la thermoformat.la \
	tinkerformat.la titleformat.la \
	turbomoleformat.la unichemformat.la \
        vaspformat.la viewmolformat.la \
	xedformat.la xtcformat.la xyzformat.la yasaraformat.la zindoformat.la
if HAVE_SHARED_POINTER
formats_LTLIBRARIES += chemdrawcdxformat.la chemkinformat.la \
                          rxnformat.la rsmiformat.la
endif
if BUILD_INCHI
formats_LTLIBRARIES += inchiformat.la
endif
if BUILD_ZLIB
formats_LTLIBRARIES += pngformat.la
endif
if HAVE_REGEX
formats_LTLIBRARIES += gamessukformat.la
endif

APIInterface_la_SOURCES = APIInterface.cpp
APIInterface_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la

CSRformat_la_SOURCES = CSRformat.cpp 
CSRformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la

MCDLformat_la_SOURCES = MCDLformat.cpp 
MCDLformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la

PQSformat_la_SOURCES = PQSformat.cpp 
PQSformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

acrformat_la_SOURCES = acrformat.cpp 
acrformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

adfformat_la_SOURCES = adfformat.cpp ArrayList.c #KReader.c
adfformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

alchemyformat_la_SOURCES = alchemyformat.cpp 
alchemyformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

amberformat_la_SOURCES = amberformat.cpp 
amberformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

balstformat_la_SOURCES = balstformat.cpp 
balstformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

bgfformat_la_SOURCES = bgfformat.cpp 
bgfformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

boxformat_la_SOURCES = boxformat.cpp 
boxformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

cacaoformat_la_SOURCES = cacaoformat.cpp 
cacaoformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

cacheformat_la_SOURCES = cacheformat.cpp 
cacheformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

carformat_la_SOURCES = carformat.cpp 
carformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

cccformat_la_SOURCES = cccformat.cpp 
cccformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

chem3dformat_la_SOURCES = chem3dformat.cpp 
chem3dformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

chemdrawctformat_la_SOURCES = chemdrawct.cpp 
chemdrawctformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

if HAVE_SHARED_POINTER
chemdrawcdxformat_la_SOURCES = chemdrawcdx.cpp 
chemdrawcdxformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 
endif

chemtoolformat_la_SOURCES = chemtoolformat.cpp 
chemtoolformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

if HAVE_SHARED_POINTER
chemkinformat_la_SOURCES = chemkinformat.cpp 
chemkinformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 
endif

cifformat_la_SOURCES = cifformat.cpp 
cifformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

copyformat_la_SOURCES = copyformat.cpp 
copyformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

crkformat_la_SOURCES = crkformat.cpp 
crkformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

cssrformat_la_SOURCES = cssrformat.cpp 
cssrformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

dmolformat_la_SOURCES = dmolformat.cpp 
dmolformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

fastsearchformat_la_SOURCES = fastsearchformat.cpp
fastsearchformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

fastaformat_la_SOURCES = fastaformat.cpp 
fastaformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

fchkformat_la_SOURCES = fchkformat.cpp 
fchkformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

featformat_la_SOURCES = featformat.cpp 
featformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

fhformat_la_SOURCES = fhformat.cpp 
fhformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la

fingerprintformat_la_SOURCES = fingerprintformat.cpp 
fingerprintformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la

freefracformat_la_SOURCES = freefracformat.cpp 
freefracformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

gamessformat_la_SOURCES = gamessformat.cpp 
gamessformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

if HAVE_REGEX
gamessukformat_la_SOURCES = gamessukformat.cpp 
gamessukformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 
endif

gaussformat_la_SOURCES = gaussformat.cpp 
gaussformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

gausscubeformat_la_SOURCES = gausscubeformat.cpp 
gausscubeformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

gausszmatformat_la_SOURCES = gausszmatformat.cpp 
gausszmatformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la

ghemicalformat_la_SOURCES = ghemicalformat.cpp 
ghemicalformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

gromos96format_la_SOURCES = gromos96format.cpp 
gromos96format_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

hinformat_la_SOURCES = hinformat.cpp 
hinformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

inchiformat_la_SOURCES = inchiformat.cpp  getinchi.cpp
if BUILD_INCHI
inchiformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la inchi102/libinchi.la
else
inchiformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la
inchiformat_la_LIBADD = -linchi
endif

jaguarformat_la_SOURCES = jaguarformat.cpp 
jaguarformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

mdlformat_la_SOURCES = mdlformat.cpp  ../alias.cpp
mdlformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

mmcifformat_la_SOURCES = mmcifformat.cpp 
mmcifformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

mmodformat_la_SOURCES = mmodformat.cpp 
mmodformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

mol2format_la_SOURCES = mol2format.cpp 
mol2format_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

molreportformat_la_SOURCES = molreport.cpp 
molreportformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

moldenformat_la_SOURCES = moldenformat.cpp 
moldenformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

mopacformat_la_SOURCES = mopacformat.cpp 
mopacformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

molproformat_la_SOURCES = molproformat.cpp 
molproformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

mpdformat_la_SOURCES = mpdformat.cpp 
mpdformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

mpqcformat_la_SOURCES = mpqcformat.cpp 
mpqcformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

msiformat_la_SOURCES = msiformat.cpp 
msiformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

msmsformat_la_SOURCES = msmsformat.cpp 
msmsformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

nwchemformat_la_SOURCES = nwchemformat.cpp 
nwchemformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

pcmodelformat_la_SOURCES = pcmodelformat.cpp 
pcmodelformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

pdbformat_la_SOURCES = pdbformat.cpp 
pdbformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

povrayformat_la_SOURCES = povrayformat.cpp 
povrayformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

pqrformat_la_SOURCES = pqrformat.cpp 
pqrformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

opendxformat_la_SOURCES = opendxformat.cpp
opendxformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

outformat_la_SOURCES = outformat.cpp
outformat_la_LDFLAGS = -module -avoid-version -no-undefined  $(win32ldflags) ../libopenbabel.la

qchemformat_la_SOURCES = qchemformat.cpp 
qchemformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32flags) ../libopenbabel.la

if BUILD_ZLIB
pngformat_la_SOURCES = pngformat.cpp 
pngformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 
pngformat_la_LIBADD = -lz
endif

reportformat_la_SOURCES = reportformat.cpp 
reportformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

if HAVE_SHARED_POINTER
rsmiformat_la_SOURCES = rsmiformat.cpp
rsmiformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la

rxnformat_la_SOURCES = rxnformat.cpp
rxnformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 
endif

shelxformat_la_SOURCES = shelxformat.cpp 
shelxformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

smilesformat_la_SOURCES = smilesformat.cpp 
smilesformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

thermoformat_la_SOURCES = thermoformat.cpp 
thermoformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

tinkerformat_la_SOURCES = tinkerformat.cpp 
tinkerformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

titleformat_la_SOURCES = titleformat.cpp 
titleformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

turbomoleformat_la_SOURCES = turbomoleformat.cpp 
turbomoleformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la

unichemformat_la_SOURCES = unichemformat.cpp 
unichemformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

vaspformat_la_SOURCES = vaspformat.cpp
vaspformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la

viewmolformat_la_SOURCES = viewmolformat.cpp 
viewmolformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

xedformat_la_SOURCES = xedformat.cpp 
xedformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

xtcformat_la_SOURCES = xtcformat.cpp 
xtcformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

xyzformat_la_SOURCES = xyzformat.cpp 
xyzformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

yasaraformat_la_SOURCES = yasaraformat.cpp 
yasaraformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

zindoformat_la_SOURCES = zindoformat.cpp 
zindoformat_la_LDFLAGS = -module -avoid-version -no-undefined $(win32ldflags) ../libopenbabel.la 

endif
