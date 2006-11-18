import os
import unittest
import pybel

class Test_readstring(unittest.TestCase):
    """Test the ability to read and write to a string"""
    def setUp(self):
        self.mol = pybel.readstring("smi", "CCCC")
        
    def accesstest(self):
        # Should raise AttributeError
        return self.mol.nosuchname

    def testformaterror(self):
        """Test that invalid formats raise an error"""
        self.assertRaises(ValueError, pybel.readstring, "noel", "jkjk")
    
    def testgetprops(self):
        """Get the values of the properties."""
        test = { 'dim':0, 'spin':1, 'energy': 0.0,
                 'charge':0, 'flags':514, 'formula': 'C4H10',
                 'mod':0 }
        result = {}
        for attr in self.mol._getmethods:
            result[attr] = getattr(self.mol, attr)
            if attr in test:
                assert result[attr] == test[attr]
        assert abs(result['exactmass']-58.078) < 0.001
        assert abs(result['molwt']-58.121) < 0.003
        self.assertEqual(len(self.mol.atoms), 4)
        self.assertRaises(AttributeError, self.accesstest)

    def testconversion(self):
        """Convert to mol2"""
        as_mol2 = self.mol.write("mol2")
        test = """@<TRIPOS>MOLECULE
*****
 4 3 0 0 0
SMALL
GASTEIGER
Energy = 0

@<TRIPOS>ATOM
      1 C           0.0000    0.0000    0.0000 C.3     1  LIG1        0.0000
      2 C           0.0000    0.0000    0.0000 C.3     1  LIG1        0.0000
      3 C           0.0000    0.0000    0.0000 C.3     1  LIG1        0.0000
      4 C           0.0000    0.0000    0.0000 C.3     1  LIG1        0.0000
@<TRIPOS>BOND
     1     1     2    1
     2     2     3    1
     3     3     4    1

"""
        self.assertEqual(as_mol2, test)

    def teststringrepr(self):
        """Test the string representation of a molecule"""
        self.assertEqual(str(self.mol).strip(), "CCCC")

class Test_readfile(unittest.TestCase):
    """Test the ability to read and write to a file"""
    def setUp(self):
        self.mols = [mol for mol in pybel.readfile("sdf", "head.sdf")]

    def testread(self):
        """Is the right number of molecules read from the file?"""
        self.assertEqual(len(self.mols), 2)

    def formaterror(self):
        mol = pybel.readfile("noel", "head.sdf").next()
    
    def testformaterror(self):
        """Test that invalid formats raise an error"""
        self.assertRaises(ValueError, self.formaterror)
    
    def testconversion(self):
        """Convert to smiles"""
        as_smi = [mol.write("smi").split("\t")[0] for mol in self.mols]
        test = ['O=C1C=CC(=O)C=C1C', 'c1cccc2c1nc(SSc1nc3ccccc3s1)s2']
        self.assertEqual(as_smi, test)

    def test_singletofile(self):
        """Test the molecule.write() method"""
        mol = self.mols[0]
        mol.write("smi", "testoutput.txt")
        test = ['O=C1C=CC(=O)C=C1C\tNSC 1\n']
        filecontents = open("testoutput.txt", "r").readlines()
        self.assertEqual(filecontents, test)
        self.assertRaises(IOError, mol.write, "smi", "testoutput.txt")
        os.remove("testoutput.txt")
        self.assertRaises(ValueError, mol.write, "noel", "testoutput.txt")
    
    def test_multipletofile(self):
        """Test the Outputfile class"""
        self.assertRaises(ValueError, pybel.Outputfile, "noel", "testoutput.txt")
        outputfile = pybel.Outputfile("smi", "testoutput.txt")
        for mol in self.mols:
            outputfile.write(mol)
        self.assertRaises(IOError, pybel.Outputfile, "smi", "testoutput.txt")
        filecontents = open("testoutput.txt", "r").readlines()
        os.remove("testoutput.txt")
        test = ['O=C1C=CC(=O)C=C1C\tNSC 1\n', 'c1cccc2c1nc(SSc1nc3ccccc3s1)s2\tNSC 2\n']
        self.assertEqual(filecontents, test)

class Test_atoms(unittest.TestCase):
    """Testing some of the atom code"""
    def setUp(self):
        self.mol = pybel.readfile("sdf", "head.sdf").next()
        self.atom = self.mol.atoms[0]

    def testiteration(self):
        """Test the ability to iterate over the atoms"""
        atoms = [atom for atom in self.mol]
        self.assertEqual(len(atoms), 15)

    def accesstest(self):
        # Should raise AttributeError
        return self.atom.nosuchname

    def testattributes(self):
        """Get the values of some properties"""
        self.assertRaises(AttributeError, self.accesstest)
        self.assert_(abs(self.atom.coords[0]-0.0021) < 0.0001)

    def teststringrepr(self):
        """Test the string representation of the Atom"""
        test = "Atom: 8 (0.0020999999999999999, -0.0041000000000000003, 0.002)"
        self.assertEqual(str(self.atom), test)

class Test_smarts(unittest.TestCase):
    """Test the Smarts object"""
    def setUp(self):
        self.mol = pybel.readstring("smi", "CCN(CC)CC")

    def testmatching(self):
        """Searching for ethyl groups in triethylamine"""
        smarts = pybel.Smarts("[#6][#6]")
        ans = smarts.findall(self.mol)
        self.assertEqual(len(ans), 3)
    
class Test_cornercases(unittest.TestCase):
    """Test some corner cases"""
    def testemptymol(self):
        """Test the creation of an empty Molecule"""
        mol = pybel.Molecule()
        self.assertEqual(mol.molwt, 0)
        self.assertEqual(len(mol.atoms), 0)
    
    def testemptyatom(self):
        """Test the creation of an empty Atom"""
        atom = pybel.Atom()
        self.assertEqual(atom.atomicnum, 0)
        
if __name__=="__main__":
    testgroups = [Test_readstring, Test_readfile, Test_cornercases, Test_atoms, Test_smarts]
    for testgroup in testgroups:
        print "\n=== %s ===" % testgroup.__doc__
        suite = unittest.makeSuite(testgroup)
        unittest.TextTestRunner(verbosity=2).run(suite)
