import os
import sys
import unittest

here = os.path.split(__file__)[0]

## In Python 3.x, generators have a __next__() method
## instead of a next() method
nextmethod = "next"
ispy2 = True
if sys.version_info[0] >= 3:
    nextmethod = "__next__"
    ispy2 = False
try:
    test = os.write
    try:
        from cinfony import pybel, rdkit, cdk
    except ImportError:
        cinfony = None
    try:
        from openbabel import pybel
        rdkit = cdk = None
    except ImportError:
        pybel = None
except AttributeError:
    from cinfony import cdk
    pybel = rdkit = None
try:
    set
except NameError:
    from sets import Set as set

class myTestCase(unittest.TestCase):
    """Additional methods not present in Jython 2.2"""
    # Taken from unittest.py in Python 2.5 distribution
    def assertFalse(self, expr, msg=None):
        "Fail the test if the expression is true."
        if expr: raise self.failureException(msg)
    def assertTrue(self, expr, msg=None):
        """Fail the test unless the expression is true."""
        if not expr: raise self.failureException(msg)
    def assertAlmostEqual(self, first, second, places=7, msg=None):
        """Fail if the two objects are unequal as determined by their
           difference rounded to the given number of decimal places
           (default 7) and comparing to zero.

           Note that decimal places (from zero) are usually not the same
           as significant digits (measured from the most signficant digit).
        """
        if round(second-first, places) != 0:
            raise self.failureException(
                  (msg or '%r != %r within %r places' % (first, second, places)))

class TestToolkit(myTestCase):

    def setUp(self):
        self.mols = [self.toolkit.readstring("smi", "CCCC"),
                     self.toolkit.readstring("smi", "CCCN")]
        self.head = list(self.toolkit.readfile("sdf", os.path.join(here, "head.sdf")))
        self.atom = self.head[0].atoms[1]

    def testattributes(self):
        """Test attributes like informats, descs and so on"""
        informats, outformats = self.toolkit.informats, self.toolkit.outformats
        self.assertNotEqual(len(list(self.toolkit.informats.keys())), 0)
        self.assertNotEqual(len(list(self.toolkit.outformats.keys())), 0)
        self.assertNotEqual(len(self.toolkit.descs), 0)
        self.assertNotEqual(len(self.toolkit.forcefields), 0)
        self.assertNotEqual(len(self.toolkit.fps), 0)

    def FPaccesstest(self):
        # Should raise AttributeError
        return self.mols[0].calcfp().nosuchname

    def testFPTanimoto(self):
        """Test the calculation of the Tanimoto coefficient"""
        fps = [x.calcfp() for x in self.mols]
        self.assertEqual(fps[0] | fps[1], self.tanimotoresult)

    def testFPstringrepr(self):
        """Test the string representation and corner cases."""
        self.assertRaises(ValueError, self.mols[0].calcfp, "Nosuchname")
        self.assertRaises(AttributeError, self.FPaccesstest)
        r = str(self.mols[0].calcfp())
        t = r.split(", ")
        self.assertEqual(len(t), self.Nfpbits)

    def testFPbits(self):
        """Test whether the bits are set correctly."""
        bits = [x.calcfp().bits for x in self.mols]
        self.assertEqual(len(bits[0]), self.Nbits)
        bits = [set(x) for x in bits]
        # Calculate the Tanimoto coefficient the old-fashioned way
        tanimoto = len(bits[0] & bits[1]) / float(len(bits[0] | bits[1]))
        self.assertEqual(tanimoto, self.tanimotoresult)

    def RSaccesstest(self):
        # Should raise AttributeError
        return self.mols[0].nosuchname

    def testRSformaterror(self):
        """Test that invalid formats raise an error"""
        self.assertRaises(ValueError, self.toolkit.readstring, "noel", "jkjk")
        self.assertRaises(IOError, self.toolkit.readstring, "smi", "&*)(%)($)")

    def testselfconversion(self):
        """Test that the toolkit can eat its own dog-food."""
        newmol = self.toolkit.Molecule(self.head[0])
        self.assertEqual(newmol._exchange,
                         self.head[0]._exchange)
        newmol = self.toolkit.Molecule(self.mols[0])
        self.assertEqual(newmol._exchange,
                         self.mols[0]._exchange)

    def testLocalOpt(self):
        """Test that local optimisation affects the coordinates"""
        oldcoords = self.head[0].atoms[0].coords
        self.head[0].localopt()
        newcoords = self.head[0].atoms[0].coords
        self.assertNotEqual(oldcoords, newcoords)
        # Make sure that make3D() is called for molecules without coordinates
        mol = self.mols[0]
        mol.localopt()
        self.assertNotEqual(mol.atoms[3].coords, (0., 0., 0.))

    def testMake2D(self):
        """Test that 2D coordinate generation does something"""
        mol = self.mols[1]
        mol.make2D()
        self.assertNotEqual(mol.atoms[2].coords, (0., 0., 0.))
        self.assertEqual(mol.atoms[2].coords[2], 0.)

    def testMake3D(self):
        """Test that 3D coordinate generation does something"""
        mol = self.mols[0]
        mol.make3D()
        self.assertNotEqual(mol.atoms[3].coords, (0., 0., 0.))

    def testDraw(self):
        """Create a 2D depiction"""
        self.mols[0].draw(show=False,
                          filename="%s.png" % self.toolkit.__name__)
        self.mols[0].draw(show=False) # Just making sure that it doesn't raise an Error
        self.mols[0].draw(show=False, update=True)
        coords = [x.coords for x in self.mols[0].atoms[0:2]]
        self.assertNotEqual(coords, [(0., 0., 0.), (0., 0., 0.)])
        self.mols[0].draw(show=False, usecoords=True,
                          filename="%s_b.png" % self.toolkit.__name__)

    def testRSgetprops(self):
        """Get the values of the properties."""
        # self.assertAlmostEqual(self.mols[0].exactmass, 58.078, 3)
        # Only OpenBabel has a working exactmass
        # CDK doesn't include implicit Hs when calculating the molwt
        self.assertAlmostEqual(self.mols[0].molwt, 58.12, 2)
        self.assertEqual(len(self.mols[0].atoms), 4)
        self.assertRaises(AttributeError, self.RSaccesstest)

    def testRSconversiontoMOL(self):
        """Convert to mol"""
        as_mol = self.mols[0].write("mol")
        test = """
 OpenBabel04220815032D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0
  1  2  1  0  0  0
  2  3  1  0  0  0
  3  4  1  0  0  0
M  END
"""
        data, result = test.split("\n"), as_mol.split("\n")
        self.assertEqual(len(data), len(result))
        self.assertEqual(data[-2], result[-2].rstrip()) # M  END

    def testRSstringrepr(self):
        """Test the string representation of a molecule"""
        self.assertEqual(str(self.mols[0]).strip(), "CCCC")

    def testRFread(self):
        """Is the right number of molecules read from the file?"""
        self.assertEqual(len(self.mols), 2)

    def RFreaderror(self):
        mol = getattr(self.toolkit.readfile("sdf", "nosuchfile.sdf"), nextmethod)()

    def testRFmissingfile(self):
        """Test that reading from a non-existent file raises an error."""
        self.assertRaises(IOError, self.RFreaderror)

    def RFformaterror(self):
        mol = getattr(self.toolkit.readfile("noel", "head.sdf"), nextmethod)()

    def testRFformaterror(self):
        """Test that invalid formats raise an error"""
        self.assertRaises(ValueError, self.RFformaterror)

    def RFunitcellerror(self):
        unitcell = self.mols[0].unitcell

    def testRFunitcellerror(self):
        """Test that accessing the unitcell raises an error"""
        self.assertRaises(AttributeError, self.RFunitcellerror)

    def testRFconversion(self):
        """Convert to smiles"""
        as_smi = [mol.write("smi").split("\t")[0] for mol in self.mols]
        ans = []
        for smi in as_smi:
            t = list(smi)
            t.sort()
            ans.append("".join(t))
        test = ['CCCC', 'CCCN']
        self.assertEqual(ans, test)

    def testRFsingletofile(self):
        """Test the molecule.write() method"""
        mol = self.mols[0]
        mol.write("smi", "testoutput.txt")
        test = 'CCCC'
        input = open("testoutput.txt", "r")
        filecontents = input.readlines()[0].split("\t")[0].strip()
        input.close()
        self.assertEqual(filecontents, test)
        self.assertRaises(IOError, mol.write, "smi", "testoutput.txt")
        os.remove("testoutput.txt")
        self.assertRaises(ValueError, mol.write, "noel", "testoutput.txt")

    def testRFoutputfile(self):
        """Test the Outputfile class"""
        self.assertRaises(ValueError, self.toolkit.Outputfile, "noel", "testoutput.txt")
        with self.toolkit.Outputfile("sdf", "testoutput.txt") as outputfile:
            for mol in self.head:
                outputfile.write(mol)
        self.assertRaises(IOError, outputfile.write, mol)
        self.assertRaises(IOError, self.toolkit.Outputfile, "sdf", "testoutput.txt")
        input = open("testoutput.txt", "r")
        numdollar = len([x for x in input.readlines()
                         if x.rstrip() == "$$$$"])
        input.close()
        os.remove("testoutput.txt")
        self.assertEqual(numdollar, 2)

    def RFdesctest(self):
        # Should raise ValueError
        self.mols[0].calcdesc("BadDescName")

    def testRFdesc(self):
        """Test the descriptors"""
        if self.toolkit.__name__ == "cinfony.cdk":
            # For the CDK, you need to call addh()
            # or some descriptors will be incorrectly calculated
            # (even those that are supposed to be immune like TPSA)
            self.mols[1].addh()
        desc = self.mols[1].calcdesc()
        self.assertTrue(len(desc) > 3)
        self.assertAlmostEqual(desc[self.tpsaname], 26.02, 2)
        self.assertRaises(ValueError, self.RFdesctest)

    def MDaccesstest(self):
        # Should raise KeyError
        return self.head[0].data['noel']

    def testMDaccess(self):
        """Change the value of a field"""
        data = self.head[0].data
        self.assertRaises(KeyError, self.MDaccesstest)
        data['noel'] = 'testvalue'
        self.assertEqual(data['noel'], 'testvalue')
        newvalues = {'hey':'there', 'yo':1}
        data.update(newvalues)
        self.assertEqual(data['yo'], '1')
        self.assertTrue('there' in data.values())

    def testMDglobalaccess(self):
        """Check out the keys"""
        data = self.head[0].data
        self.assertFalse('Noel' in data)
        self.assertEqual(len(data), len(self.datakeys))
        for key in data:
            self.assertEqual(key in self.datakeys, True)
        r = repr(data)
        self.assertTrue(r[0]=="{" and r[-2:]=="'}", r)

    def testMDdelete(self):
        """Delete some keys"""
        data = self.head[0].data
        self.assertTrue('NSC' in data)
        del data['NSC']
        self.assertFalse('NSC' in data)
        data.clear()
        self.assertEqual(len(data), 0)

    def testAiteration(self):
        """Test the ability to iterate over the atoms"""
        atoms = [atom for atom in self.head[0]]
        self.assertEqual(len(atoms), self.Natoms)

    def Atomaccesstest(self):
        # Should raise AttributeError
        return self.atom.nosuchname

    def testAattributes(self):
        """Get the values of some properties"""
        self.assertRaises(AttributeError, self.Atomaccesstest)
        self.assertAlmostEqual(self.atom.coords[0], -0.0691, 4)

    def testAstringrepr(self):
        """Test the string representation of the Atom"""
        test = "Atom: 8 (-0.07 5.24 0.03)"
        self.assertEqual(str(self.atom), test)

    def invalidSMARTStest(self):
        # Should raise IOError
        return self.toolkit.Smarts("[#NOEL][#NOEL]")

    def testSMARTS(self):
        """Searching for ethyl groups in triethylamine"""
        mol = self.toolkit.readstring("smi", "CCN(CC)CC")
        smarts = self.toolkit.Smarts("[#6][#6]")
        ans = smarts.findall(mol)
        self.assertEqual(len(ans), 3)
        self.toolkit.ob.obErrorLog.SetOutputLevel(self.toolkit.ob.obError)
        self.assertRaises(IOError, self.invalidSMARTStest)
        self.toolkit.ob.obErrorLog.SetOutputLevel(self.toolkit.ob.obWarning)

    def testAddh(self):
        """Adding and removing hydrogens"""
        self.assertEqual(len(self.mols[0].atoms),4)
        self.mols[0].addh()
        self.assertEqual(len(self.mols[0].atoms),14)
        self.mols[0].removeh()
        self.assertEqual(len(self.mols[0].atoms),4)

class TestPybel(TestToolkit):
    toolkit = pybel
    tanimotoresult = 1/3.
    Natoms = 15
    tpsaname = "TPSA"
    Nbits = 3
    Nfpbits = 32
    datakeys = ['NSC', 'Comment', 'OpenBabel Symmetry Classes', 'MOL Chiral Flag']

    def testFP_FP3(self):
        "Checking the results from FP3"
        fps = [x.calcfp("FP3") for x in self.mols]
        self.assertEqual(fps[0] | fps[1], 0.)

    def testunitcell(self):
        """Testing unit cell access"""
        mol = getattr(self.toolkit.readfile("cif", os.path.join(here, "hashizume.cif")), nextmethod)()
        cell = mol.unitcell
        self.assertAlmostEqual(cell.GetAlpha(), 93.0, 1)

    def testMDcomment(self):
        """Mess about with the comment field"""
        data = self.head[0].data
        self.assertEqual('Comment' in data, True)
        self.assertEqual(data['Comment'], 'CORINA 2.61 0041  25.10.2001')
        data['Comment'] = 'New comment'
        self.assertEqual(data['Comment'], 'New comment')

    def importtest(self):
        self.mols[0].draw(show=True, usecoords=True)

    def testRSconversiontoMOL2(self):
        """Convert to mol2"""
        as_mol2 = self.mols[0].write("mol2")
        test = """@<TRIPOS>MOLECULE
*****
 4 3 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C           0.0000    0.0000    0.0000 C.3     1  UNL1        0.0000
      2 C           0.0000    0.0000    0.0000 C.3     1  UNL1        0.0000
      3 C           0.0000    0.0000    0.0000 C.3     1  UNL1        0.0000
      4 C           0.0000    0.0000    0.0000 C.3     1  UNL1        0.0000
@<TRIPOS>BOND
     1     1     2    1
     2     2     3    1
     3     3     4    1
"""
        self.assertEqual(as_mol2, test)

    def testRSgetprops(self):
        """Get the values of the properties."""
        self.assertAlmostEqual(self.mols[0].exactmass, 58.078, 3)
        self.assertAlmostEqual(self.mols[0].molwt, 58.122, 3)
        self.assertEqual(len(self.mols[0].atoms), 4)
        self.assertRaises(AttributeError, self.RSaccesstest)

    def testIterators(self):
        """Check out the OB iterators"""
        numatoms = len(list(self.toolkit.ob.OBMolAtomIter(self.mols[0].OBMol)))
        self.assertEqual(numatoms, 4)

class TestOBPybelNoDraw(TestPybel):
    def testDraw(self):
       """No drawing done"""
       pass

class TestPybelWithDraw(TestPybel):

    def testDrawdependencies(self):
        "Testing the draw dependencies"
        t = self.toolkit.tk
        self.toolkit.tk = None
        self.mols[0].draw(show=False, usecoords=True,
                          filename="%s_b.png" % self.toolkit.__name__)
        self.assertRaises(ImportError,
                          self.importtest)
        self.toolkit.tk = t

        t = self.toolkit.oasa
        self.toolkit.oasa = None
        self.assertRaises(ImportError,
                          self.importtest)


class TestRDKit(TestToolkit):
    toolkit = rdkit
    tanimotoresult = 1/3.
    Natoms = 9
    tpsaname = "TPSA"
    Nbits = 12
    Nfpbits = 64
    datakeys = ['NSC']


class TestCDK(TestToolkit):
    toolkit = cdk
    tanimotoresult = 0.375
    Natoms = 15
    tpsaname = "tpsa"
    Nbits = 4
    Nfpbits = 4 # The CDK uses a true java.util.Bitset
    datakeys = ['NSC', 'Remark', 'Title']

    def testSMARTS(self):
        """No SMARTS testing done"""
        pass

    def testLocalOpt(self):
        """No local opt testing done"""
        pass

    def testMake2D(self):
        """No 2D coordinate generation done"""
        pass

    def testMake3D(self):
        """No 3D coordinate generation done"""
        pass

    def testRSgetprops(self):
        """Get the values of the properties."""
        # self.assertAlmostEqual(self.mols[0].exactmass, 58.078, 3)
        # Only OpenBabel has a working exactmass
        # CDK doesn't include implicit Hs when calculating the molwt
        self.mols[0].addh()
        self.assertAlmostEqual(self.mols[0].molwt, 58.12, 2)
        self.assertEqual(len(self.mols[0].atoms), 14)
        self.assertRaises(AttributeError, self.RSaccesstest)

if __name__=="__main__":
    # Tidy up
    if os.path.isfile("testoutput.txt"):
        os.remove("testoutput.txt")

    #testcases = [TestPybel, TestCDK, TestRDKit]
    # testcases = [TestCDK]
    # testcases = [TestPybel]
    # testcases = [TestRDKit]
    testcases = [TestPybel]
    for testcase in testcases:
        sys.stdout.write("\n\n\nTESTING %s\n%s\n\n\n" % (testcase.__name__, "== "*10))
        myunittest = unittest.defaultTestLoader.loadTestsFromTestCase(testcase)
        unittest.TextTestRunner(verbosity=2).run(myunittest)
