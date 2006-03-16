import pyopenbabel as pob
import openbabel as ob

def testwrapper():
    """Test how well the wrapper is wrapped."""

    print "\n******* How well does pyopenbabel wrap openbabel? *********\n\n"

    print "*** Atom methods ***\n"

    getmethods = set([x for x in ob.OBAtom.__dict__.keys() if x.startswith("Get")])
    wrapped = set(pob.Atom._getmethods.values())

    print "(1) 'Get methods' wrapped by SWIG and included in pyopenbabel:",sorted(getmethods & wrapped)
    print "\n(2) 'Get methods' wrapped by SWIG but not in pyopenbabel:",sorted(getmethods-wrapped)
    print "\n(3) Mistakes: 'Get methods' in pyopenbabel but not in SWIG:",sorted(wrapped-getmethods)


    print "\n\n*** Molecule methods ***\n"

    getmethods = set([x for x in ob.OBMol.__dict__.keys() if x.startswith("Get")])
    wrapped = set(pob.Molecule._getmethods.values())

    print "(1) 'Get methods' wrapped by SWIG and included in pyopenbabel:",sorted(getmethods & wrapped)
    print "\n(2) 'Get methods' wrapped by SWIG but not in pyopenbabel:",sorted(getmethods-wrapped)
    print "\n(3) Mistakes: 'Get methods' in pyopenbabel but not in SWIG:",sorted(wrapped-getmethods)

def testsd(sdfile):
    """Read in an sdfile."""
    sd = open(sdfile,"r").read()
    mol = pob.readstring("sd",sd)
    for x,y in mol.__dict__.iteritems():
        print x,y

if __name__=="__main__":
    # testwrapper()
    testsd("3d.head.sdf")
