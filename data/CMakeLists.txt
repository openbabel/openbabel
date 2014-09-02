# Open Babel data files

set(to_install
    aromatic.txt
    atomization-energies.txt
    atomtyp.txt
    babel_povray3.inc
    bondtyp.txt
    element.txt
    eqeqIonizations.txt
    fragments.txt
    gaff.dat
    gaff.prm
    ghemical.prm
    isotope-small.txt
    isotope.txt
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
    ringtyp.txt
    SMARTS_InteLigand.txt
    space-groups.txt
    superatom.txt
    svgformat.script
    templates.sdf
    torlib.txt
    types.txt
    UFF.prm
)

if(NOT MSVC)
  install(FILES ${to_install} DESTINATION share/openbabel/${BABEL_VERSION})
else(NOT MSVC)
  install(FILES ${to_install} DESTINATION bin/data)
endif(NOT MSVC)
