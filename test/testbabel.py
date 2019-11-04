"""Test OpenBabel executables from Python

Note: Python bindings not used

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pytest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../test/testbabel.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import os
import re
import sys
import unittest

from subprocess import CalledProcessError, Popen, PIPE, check_output, STDOUT

def run_exec(*args):
    """Run one of OpenBabel's executables

    With two arguments (stdin, commandline) it pipes
    the stdin through the executable.

    Example: run_exec("CC(=O)Cl", "obabel -ismi -oinchi")

    Return a tuple (stdout, stderr)
    """
    if len(args) == 2:
        text, commandline = args
    elif len(args) == 1:
        text, commandline = "", args[0]
    else:
        raise Exception("One or two arguments expected")

    if sys.platform.startswith("win"):
        broken = commandline.split()
        exe = executable(broken[0])
    else:
        broken = commandline.encode().split()
        exe = executable(broken[0].decode())
    # Note that bufsize = -1 means default buffering
    # Without this, it's unbuffered and it takes 10x longer on MacOSX
    if text:
        p = Popen([exe] + broken[1:],
                  stdin=PIPE, stdout=PIPE, stderr=PIPE, bufsize=1)
        stdout, stderr = p.communicate(text.encode())
    else:
        p = Popen([exe] + broken[1:],
                  stdout=PIPE, stderr=PIPE, bufsize=-1)
        stdout, stderr = p.communicate()

    return stdout.decode(), stderr.decode()

def executable(name):
    """Return the full path to an executable"""
    suffix = ""
    folder = "bin"
    if sys.platform == "win32":
        suffix = ".exe"
        folder = os.path.join(folder, "Release")
    return os.path.join("..", folder, name + suffix)

def log(text):
    """Convenience function for debugging tests

    The log file (log.txt) is created in build/test
    """
    with open("log.txt", "a") as output:
        output.write(text + "\n")

class BaseTest(unittest.TestCase):
    """A base class for test classes that adds additional
    test methods"""

    def canFindExecutable(self, name):
        fullpath = executable(name)
        self.assertTrue(os.path.isfile(fullpath),
                        "'%s' executable not found at %s" % (name, fullpath))

    def canFindFile(self, filename):
        self.assertTrue(os.path.isfile(filename),
                        "Cannot find the file '%s'" % filename)

    def getTestFile(self, filename):
        here = sys.path[0]
        fullpath = os.path.join(here, "files", filename)
        self.canFindFile(fullpath)
        return fullpath

    def assertConverted(self, stderr, N):
        """Assert that N molecules were converted."""
        pat = "(-?\d.*) molecule(?:s?) converted"
        lines = stderr.split("\n")
        convertedline =  [line for line in lines if re.match(pat, line)]
        if len(convertedline) == 0:
            self.fail("Cannot find the number of molecules converted")
        conversion_no = int(re.findall(pat, convertedline[0])[0])
        self.assertEqual(N, conversion_no,
                         "Number of molecules converted is %d "
                         "but should be %d" % (conversion_no, N))

class TestOBabel(BaseTest):
    """A series of tests relating to the obabel executable"""

    def testNoInput(self):
        """Ensure that this does not segfault (PR#1818)"""
        self.canFindExecutable("obabel")
        output, error = run_exec("obabel -i")

    def testSMItoInChI(self):
        self.canFindExecutable("obabel")
        output, error = run_exec("CC(=O)Cl", "obabel -ismi -oinchi")
        self.assertEqual(output.rstrip(), "InChI=1S/C2H3ClO/c1-2(3)4/h1H3")

    def testRSMItoRSMI(self):
        # Check possible combinations of missing rxn components
        data = ["O>N>S", "O>>S", "O>N>", "O>>",
                ">N>S", ">>S", ">>"]
        for rsmi in data:
            output, error = run_exec('obabel -:%s -irsmi -orsmi' % rsmi)
            self.assertEqual(output.rstrip(), rsmi)
        # Check handling of invalid rxn components
        data = ["Noel>N>S", "O>Noel>S", "O>N>Noel"]
        errors = ["reactant", "agent", "product"]
        for rsmi, error in zip(data, errors):
            output, errormsg = run_exec('obabel -:%s -irsmi -orsmi' % rsmi)
            self.assertTrue(error in errormsg)

    def sort(self, rsmi):
        # TODO: Change OBMol.Separate to preserve the order. This
        # function shouldn't be necessary.
        tmp = rsmi.split(">")
        r = ".".join(sorted(tmp[0].split(".")))
        p = ".".join(sorted(tmp[2].split(".")))
        return "%s>%s>%s" % (r, tmp[1], p)

    def testRingClosures(self):
        # Test positives
        data = ["c1ccccc1", "c%11ccccc%11", "c%(1)ccccc%(1)", "c%(51)ccccc%51",
                "c%(99999)ccccc%(99999)"]
        for smi in data:
            output, error = run_exec("obabel -:%s -osmi" % smi)
            self.assertEqual("c1ccccc1", output.rstrip())
        # Test negatives
        data = ["c%1ccccc%1", "c%a1cccc%a1", "c%(a1)ccccc%(a1)",
                "c%(000001)ccccc%(000001)", "c%(51)ccccc%(15)"]
        for smi in data:
            output, error = run_exec("obabel -:%s -osmi" % smi)
            self.assertTrue("0 molecules converted" in error)
        # Now test writing of %(NNN) notation
        output, error = run_exec("obabel %s -osmi" % self.getTestFile("102Uridine.smi"))
        self.assertTrue("%(100)" in output)

    def testPDBQT(self):
        self.canFindExecutable("obabel")
        pdb = '''ATOM     77  N   TYR A   5      35.078  50.693  67.193  1.00  0.00           N  
ATOM     78  CA  TYR A   5      35.195  51.589  66.041  1.00  0.00           C  
ATOM     79  C   TYR A   5      33.792  51.581  65.423  1.00  0.00           C  
ATOM     80  O   TYR A   5      33.362  50.580  64.852  1.00  0.00           O  
ATOM     81  CB  TYR A   5      36.233  51.055  65.045  1.00  0.00           C  
ATOM     82  CG  TYR A   5      36.534  51.990  63.891  1.00  0.00           C  
ATOM     83  CD1 TYR A   5      37.369  53.085  64.053  1.00  0.00           C  
ATOM     84  CD2 TYR A   5      35.949  51.786  62.644  1.00  0.00           C  
ATOM     85  CE1 TYR A   5      37.618  53.971  62.993  1.00  0.00           C  
ATOM     86  CE2 TYR A   5      36.180  52.653  61.580  1.00  0.00           C  
ATOM     87  CZ  TYR A   5      37.010  53.742  61.763  1.00  0.00           C  
ATOM     88  OH  TYR A   5      37.199  54.607  60.730  1.00  0.00           O  
ATOM     89  H   TYR A   5      35.302  49.703  67.057  1.00  0.00           H  
ATOM     90  HA  TYR A   5      35.457  52.597  66.383  1.00  0.00           H  
ATOM     91  HB2 TYR A   5      35.921  50.094  64.629  1.00  0.00           H  
ATOM     92  HB3 TYR A   5      37.173  50.847  65.572  1.00  0.00           H  
ATOM     93  HD1 TYR A   5      37.853  53.262  65.012  1.00  0.00           H  
ATOM     94  HD2 TYR A   5      35.274  50.955  62.515  1.00  0.00           H  
ATOM     95  HE1 TYR A   5      38.313  54.794  63.081  1.00  0.00           H  
ATOM     96  HE2 TYR A   5      35.711  52.433  60.627  1.00  0.00           H  
ATOM     97  HH  TYR A   5      36.875  54.171  59.926  1.00  0.00           H  
END
'''
        pdbqt = '''REMARK  Name = 
REMARK  5 active torsions:
REMARK  status: ('A' for Active; 'I' for Inactive)
REMARK    1  A    between atoms: N_1  and  CA_2
REMARK    2  A    between atoms: CA_2  and  C_3
REMARK    3  A    between atoms: CA_2  and  CB_5
REMARK    4  A    between atoms: CB_5  and  CG_6
REMARK    5  A    between atoms: CZ_11  and  OH_12
REMARK                            x       y       z     vdW  Elec       q    Type
REMARK                         _______ _______ _______ _____ _____    ______ ____
ROOT
ATOM      1  CG  TYR A   5      36.534  51.990  63.891  0.00  0.00    +0.000 A 
ATOM      2  CD1 TYR A   5      37.369  53.085  64.053  0.00  0.00    +0.000 A 
ATOM      3  CD2 TYR A   5      35.949  51.786  62.644  0.00  0.00    +0.000 A 
ATOM      4  CE1 TYR A   5      37.618  53.971  62.993  0.00  0.00    +0.000 A 
ATOM      5  CE2 TYR A   5      36.180  52.653  61.580  0.00  0.00    +0.000 A 
ATOM      6  CZ  TYR A   5      37.010  53.742  61.763  0.00  0.00    +0.000 A 
ENDROOT
BRANCH   6   7
ATOM      7  OH  TYR A   5      37.199  54.607  60.730  0.00  0.00    +0.000 OA
ATOM      8  HH  TYR A   5      36.875  54.171  59.926  0.00  0.00    +0.000 HD
ENDBRANCH   6   7
BRANCH   1   9
ATOM      9  CB  TYR A   5      36.233  51.055  65.045  0.00  0.00    +0.000 C 
BRANCH   9  10
ATOM     10  CA  TYR A   5      35.195  51.589  66.041  0.00  0.00    +0.000 C 
BRANCH  10  11
ATOM     11  N   TYR A   5      35.078  50.693  67.193  0.00  0.00    +0.000 NA
ATOM     12  H   TYR A   5      35.302  49.703  67.057  0.00  0.00    +0.000 HD
ENDBRANCH  10  11
BRANCH  10  13
ATOM     13  C   TYR A   5      33.792  51.581  65.423  0.00  0.00    +0.000 C 
ATOM     14  O   TYR A   5      33.362  50.580  64.852  0.00  0.00    +0.000 OA
ENDBRANCH  10  13
ENDBRANCH   9  10
ENDBRANCH   1   9
TORSDOF 5
'''
        output, error = run_exec(pdb, "obabel -ipdb -opdbqt")
        self.assertEqual(output.replace("\r", ""), pdbqt.replace("\r", ""))

    def testMissingPlugins(self):
        if sys.platform.startswith("win32"):
            return
        libdir = os.environ.pop("BABEL_LIBDIR", None)
        os.environ["BABEL_LIBDIR"] = ""

        obabel = executable("obabel")
        with self.assertRaises(CalledProcessError) as cm:
            check_output('%s -:C -osmi' % obabel, shell=True, stderr=STDOUT, universal_newlines=True)
        msg = cm.exception.output
        if libdir:
            os.environ["BABEL_LIBDIR"] = libdir
        else:
            os.environ.pop("BABEL_LIBDIR")

        self.assertTrue('BABEL_LIBDIR' in msg)

    def testCOFtoCAN(self):
        self.canFindExecutable("obabel")
        listCOFnames = [
                'culgi_00',
                'culgi_01',
                'culgi_02',
                'culgi_03',
                'culgi_04',
                'culgi_05',
                'culgi_06',
                'culgi_07',
                'culgi_08',
                'culgi_09',
                'culgi_10',
                'culgi_11',
                ]
        listCANexpected = [
                '[O-]C(=O)Cc1cccc2c1oc1c(C)c(C)ccc1c2=O',
                'C=Cc1c[nH]c(=O)[nH]c1=O',
                'C[C@@H](Cc1ccc(cc1)I)[NH2+]C(C)C',
                '[NH3+]C[C@@H]1O[C@@H](Cc2c1ccc(c2O)O)[C@]12C[C@H]3C[C@@H](C2)C[C@@H](C1)C3',
                '[O-]C(=O)[C@H](CCCNC(=[NH2+])N)[NH3+]',
                'ClC(=O)C(C)(C)C',
                'CC#CCOC(=O)c1c(C)nc2c(c1N)c1CC[C@H](Cc1s2)O',
                'OC[C@H]1O[C@H](C[C@@H]1F)n1ccc(=O)[nH]c1=O',
                'OC[C@H]1O[C@H]([C@@H]([C@@H]1O)O)n1cnc2c1ncnc2N[C@H]1CCC[C@@H]1O',
                '[O-]C(=O)CC[C@H]([NH3+])C=C',
                'O=[S@@](c1nc2c([nH]1)cccc2)Cc1nccc(c1C)OCC(F)(F)F',
                'C#C[C@]1(O)CC[C@@H]2[C@]1(C)CC[C@H]1[C@H]2CCc2c1ccc(c2)O',
                ]
        for cofname, CAN in zip(listCOFnames, listCANexpected):
            coffilename = cofname + '.cof'
            if(cofname == 'culgi_06'):
                cofname = 'mol24' # Special case: 'internal name' not the same as file name
            coffile = self.getTestFile(coffilename)
            cansmi = CAN + '\t' + cofname # Expected SMILES line plus molecule name
            output, error = run_exec( "obabel -icof -ocan %s" % coffile)
            self.assertEqual(output.rstrip('\r\n'), cansmi)

    def testCOFtoMOL(self):
        self.canFindExecutable("obabel")
        listCOFnames = [
                'culgi_00',
                'culgi_01',
                'culgi_02',
                'culgi_03',
                'culgi_04',
                'culgi_05',
                'culgi_06',
                'culgi_07',
                'culgi_08',
                'culgi_09',
                'culgi_10',
                'culgi_11',
                ]
        for cofname in listCOFnames:
            coffilename = cofname + '.cof'
            if(cofname == 'culgi_06'):
                cofname = 'mol24'
            coffile = self.getTestFile(coffilename)
            output, error = run_exec( "obabel -icof -omol %s" % coffile)
            molfilename = cofname + '.mol'
            molfile = self.getTestFile(molfilename)

            # Chop up the output and the baseline files into single lines
            # Skip first two lines: first line contains Culgi version,
            # second line contains OpenBabel ID/hash
            with open(molfile, "r") as molfilehandle:
                moldata = molfilehandle.readlines()[2:]
            outdata = output.splitlines()[2:]
            self.assertEqual(len(outdata), len(moldata))
            for outline, molline in zip(outdata, moldata):
                self.assertEqual(outline.rstrip('\r\n'), molline.rstrip('\r\n'))

    def testMOLtoCOF(self):
        self.canFindExecutable("obabel")
        listMOLnames = [
                'culgi_00',
                'culgi_01',
                'culgi_02',
                'culgi_03',
                'culgi_04',
                'culgi_05',
                'mol24',
                'culgi_07',
                'culgi_08',
                'culgi_09',
                'culgi_10',
                'culgi_11',
                ]
        for molname in listMOLnames:
            molfilename = molname + '.mol'
            coffilename = molname + '_from_mol.cof'
            coffile = self.getTestFile(coffilename)
            molfile = self.getTestFile(molfilename)
            output, error = run_exec( "obabel -imol -ocof %s --partialcharge none" % molfile)

            # Chop up the output and the baseline files into single lines
            # Skip first three lines: first line contains Culgi version,
            # next two lines contain comment
            with open(coffile, "r") as coffilehandle:
                cofdata = coffilehandle.readlines()[3:]
            outdata = output.splitlines()[3:]
            self.assertEqual(len(outdata), len(cofdata))
            for outline, cofline in zip(outdata, cofdata):
                self.assertEqual(outline.rstrip('\r\n'), cofline.rstrip('\r\n'))

    def testCOFtoMOL2(self):
        self.canFindExecutable("obabel")
        listCOFnames = [
                'culgi_00',
                'culgi_01',
                'culgi_02',
                'culgi_03',
                'culgi_04',
                'culgi_05',
                'culgi_06',
                'culgi_07',
                'culgi_08',
                'culgi_09',
                'culgi_10',
                'culgi_11',
                ]
        for cofname in listCOFnames:
            coffilename = cofname + '.cof'
            if(cofname == 'culgi_06'):
                cofname = 'mol24'
            coffile = self.getTestFile(coffilename)
            output, error = run_exec( "obabel -icof -omol2 %s" % coffile)
            mol2filename = cofname + '.mol2'
            mol2file = self.getTestFile(mol2filename)

            # Chop up the output and the baseline files into single lines
            with open(mol2file, "r") as mol2filehandle:
                mol2data = mol2filehandle.readlines()
            outdata = output.splitlines()
            self.assertEqual(len(outdata), len(mol2data))
            for outline, mol2line in zip(outdata, mol2data):
                self.assertEqual(outline.rstrip('\r\n'), mol2line.rstrip('\r\n'))

    def testMOL2toCOF(self):
        self.canFindExecutable("obabel")
        listMOL2names = [
                'culgi_00',
                'culgi_01',
                'culgi_02',
                'culgi_03',
                'culgi_04',
                'culgi_05',
                'mol24',
                'culgi_07',
                'culgi_08',
                'culgi_09',
                'culgi_10',
                'culgi_11',
                ]
        for mol2name in listMOL2names:
            mol2filename = mol2name + '.mol2'
            coffilename = mol2name + '_from_mol2.cof'
            coffile = self.getTestFile(coffilename)
            mol2file = self.getTestFile(mol2filename)
            output, error = run_exec( "obabel -imol2 -ocof %s" % mol2file)

            # Chop up the output and the baseline files into single lines
            # Skip first three lines: first line contains Culgi version,
            # next two lines contain comment
            with open(coffile, "r") as coffilehandle:
                cofdata = coffilehandle.readlines()[3:]
            outdata = output.splitlines()[3:]
            self.assertEqual(len(outdata), len(cofdata))
            for outline, cofline in zip(outdata, cofdata):
                self.assertEqual(outline.rstrip('\r\n'), cofline.rstrip('\r\n'))

if __name__ == "__main__":
    unittest.main()
