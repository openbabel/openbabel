#!/bin/sh

./bin2hex.pl atomtyp.txt AtomTypeData >atomtyp.h
./bin2hex.pl bondtyp.txt BondTypeData >bondtyp.h
./bin2hex.pl phmodel.txt PhModelData PHMODELDATA >phmodeldata.h
./bin2hex.pl resdata.txt ResidueData >resdata.h
./bin2hex.pl torlib.txt TorsionDefaults >torlib.h
./bin2hex.pl types.txt TypesData >types.h
./bin2hex.pl atomization-energies.txt AtomicHeatOfFormationData > atomizationenergies.h
./bin2hex.pl space-groups.txt SpaceGroupsData > spacegroups.h
./bin2hex.pl ringtyp.txt RingTypeData > ringtyp.h
