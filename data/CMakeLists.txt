# 2008-03-15

set(to_install
    aromatic.txt
    atomtyp.txt
    bondtyp.txt
    element.txt
    fragments.txt
    ghemical.prm
    isotope-small.txt
    isotope.txt
    logp.txt
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
    mr.txt
    patterns.txt
    phmodel.txt
    plugindefines.txt
    psa.txt
    resdata.txt
    ringtyp.txt
    SMARTS_InteLigand.txt
    space-groups.txt
    space-groups.xsl # needed?
    torlib.txt
    types.txt
    UFF.prm
)

install(FILES ${to_install} DESTINATION share/openbabel/${BABEL_VERSION})
