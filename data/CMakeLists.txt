
set(to_install
    aromatic.txt bondtyp.txt isotope.txt resdata.txt
    atomtyp.txt element.txt isotope-small.txt phmodel.txt
    patterns.txt torlib.txt types.txt
    SMARTS_InteLigand.txt
    ghemical.prm
    logp.txt psa.txt mr.txt
)

install(FILES ${to_install} DESTINATION share/openbabel/${BABEL_VERSION})
