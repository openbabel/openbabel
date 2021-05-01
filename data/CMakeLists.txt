# Open Babel data files

set(to_install
    atomization-energies.txt
    atomtyp.txt
    babel_povray3.inc
    bondtyp.txt
    eem.txt
    eem2015ba.txt
    eem2015bm.txt
    eem2015bn.txt
    eem2015ha.txt
    eem2015hm.txt
    eem2015hn.txt
    eqeqIonizations.txt
    gaff.dat
    gaff.prm
    ghemical.prm
    logp.txt
    MACCS.txt
    mmff94.ff
    mmff94s.ff
    mmffang.par
    mmffbndk.par
    mmffbond.par
    mmffchg.par
    mmffdef.par
    mmffdfsb.par
    mmffoop.par
    mmffpbci.par
    mmffprop.par
    mmffstbn.par
    mmfftor.par
    mmffvdw.par
    mmffs_oop.par
    mmffs_tor.par
    mpC.txt
    mr.txt
    patterns.txt
    phmodel.txt
    plugindefines.txt
    psa.txt
    qeq.txt
    resdata.txt
    rigid-fragments.txt
    rigid-fragments-index.txt
    ring-fragments.txt
    ringtyp.txt
    SMARTS_InteLigand.txt
    space-groups.txt
    superatom.txt
    svgformat.script
    templates.sdf
    torlib.txt
    torsion.txt
    types.txt
    UFF.prm
)

if(NOT MSVC)
  install(FILES ${to_install} DESTINATION share/openbabel/${BABEL_VERSION})
else(NOT MSVC)
  install(FILES ${to_install} DESTINATION bin/data)
endif(NOT MSVC)



#input_file variable_name header
set(datainfo
"atomtyp.txt AtomTypeData atomtyp.h"
"bondtyp.txt BondTypeData bondtyp.h"
"phmodel.txt PhModelData phmodeldata.h"
"resdata.txt ResidueData resdata.h"
"torlib.txt TorsionDefaults torlib.h"
"types.txt TypesData types.h"
"atomization-energies.txt AtomicHeatOfFormationData atomizationenergies.h"
"space-groups.txt SpaceGroupsData spacegroups.h"
"ringtyp.txt RingTypeData ringtyp.h"
)

set(dataheaderfiles "")
foreach(dinfo ${datainfo})
 string(REPLACE " " ";" dlist ${dinfo})
 list(GET dlist 2 doutput)
 list(GET dlist 1 dname)
 list(GET dlist 0 dinput)

 add_custom_command(OUTPUT ${doutput}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin2hex.pl ${CMAKE_CURRENT_SOURCE_DIR}/${dinput} ${dname} ${dname} > ${doutput}
    DEPENDS ${dinput} VERBATIM)    
 set_source_files_properties(${doutput} PROPERTIES GENERATED TRUE)
 list(APPEND dataheaderfiles ${doutput})
endforeach(dinfo)

add_custom_target(data_headers DEPENDS ${dataheaderfiles})
