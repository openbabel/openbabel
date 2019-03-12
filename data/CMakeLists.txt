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
