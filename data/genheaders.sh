#!/bin/sh

./bin2hex.pl aromatic.txt AromaticData >aromatic.h
./bin2hex.pl atomtyp.txt AtomTypeData >atomtyp.h
./bin2hex.pl bondtyp.txt BondTypeData >bondtyp.h
./bin2hex.pl element.txt ElementData >element.h
./bin2hex.pl isotope-small.txt IsotopeData ISOTOPE >isotope.h
./bin2hex.pl phmodel.txt PhModelData PHMODELDATA >phmodeldata.h
./bin2hex.pl resdata.txt ResidueData >resdata.h
./bin2hex.pl torlib.txt TorsionDefaults >torlib.h
./bin2hex.pl types.txt TypesData >types.h
