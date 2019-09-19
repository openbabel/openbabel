"""Test the ability of openbabel to roundtrip various properties
(protonation, hybridization, aromaticity...) through various file formats
(sdf,pdb,mol2).

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pyroundtriptest -VV

The runtime directory is ${CMAKE_SRC_DIR}/test.

You could also "chdir" into build and run the test file directly:
python ../../test/testroundtrip.py

In this latter case, you will need to set the environment variables
PYTHONPATH, LD_LIBRARY_PATH, BABEL_LIBDIR and BABEL_DATADIR beforehand.
The CMake script does this automatically.

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import os
import sys
import unittest
import itertools
import glob
import traceback

import itertools
from multiprocessing import Pool  
from multiprocessing import TimeoutError

from contextlib import contextmanager

here = sys.path[0]
iswin = sys.platform.startswith("win")

try:
    from openbabel import openbabel as ob
except ImportError:
    ob = None

try:
    from openbabel import pybel
except ImportError:
    pybel = None
    

def molsAreSame(a, b):
  '''Return true if a and b are the same'''
  n = len(a.atoms)
  if n != len(b.atoms):
    return "Different number of atoms, %d != %d"%(n,len(b.atoms))
    
  iso = ob.OBIsomorphismMapper.GetInstance(ob.CompileMoleculeQuery(a.OBMol))
  m = ob.vpairUIntUInt() 
  iso.MapFirst(b.OBMol,m)
  
  if len(m) != n:
    badsdf = open('badsdf.sdf','wt')
    badsdf.write(a.write('sdf'))
    badsdf.write(b.write('sdf'))
    badsdf.close()
    badmol = open('badmol.mol2','wt')
    badmol.write(a.write('mol2'))
    badmol.write(b.write('mol2'))
    badmol.close()    
    return "Could only match %d of %d atoms (%d)"%(len(m),n,len(b.atoms))
  
  for (ai,bi) in m:
    aatom = a.OBMol.GetAtom(ai+1)
    batom = b.OBMol.GetAtom(bi+1)
    
    if aatom.GetAtomicNum() != batom.GetAtomicNum():
      return "Mismatched atom num, %d != %d"%(aatom.GetAtomicNum(), batom.GetAtomicNum())
    if aatom.GetHyb() != batom.GetHyb():
      return "Mismatched Hyb, %d != %d"%(aatom.GetHyb(), batom.GetHyb())
    if aatom.GetExplicitValence() != batom.GetExplicitValence():
      return "Mismatched ExplicitValence, %d != %d"%(aatom.GetExplicitValence(), batom.GetExplicitValence())
    if aatom.IsAromatic () != batom.IsAromatic ():
      return "Mismatched IsAromatic, %d != %d"%( aatom.IsAromatic(), batom.IsAromatic())
    if aatom.GetFormalCharge() != batom.GetFormalCharge():
      return "Mismatched GetFormalCharge, %d != %d"%(aatom.GetFormalCharge(), batom.GetFormalCharge())
         
  return None
  
def dumpMol(m):
    for a in m.atoms:
      A = a.OBAtom
      print('%d %d %d  %d  %d  %d  (%f,%f,%f)'%(A.GetId(), A.GetIdx(), A.GetAtomicNum(),\
        A.GetHyb(),A.GetExplicitValence(),A.IsAromatic(), \
        A.GetX(),A.GetY(),A.GetZ()))
      
def dumpBoth(a, b):
  '''Print out both molecules in tandem'''
  print("Idx  ANum  Hyb  Val  Aro    Coords")
  for (a,b) in zip(a.atoms,b.atoms):
    A = a.OBAtom
    B = b.OBAtom
    print('%d:%d  %d:%d  %d:%d  %d:%d  %d:%d  (%f,%f,%f):(%f,%f,%f)'%(A.GetIdx(),B.GetIdx(),A.GetAtomicNum(),B.GetAtomicNum(),\
      A.GetHyb(),B.GetHyb(),A.GetExplicitValence(),B.GetExplicitValence(),A.IsAromatic(),B.IsAromatic(), \
      A.GetX(),A.GetY(),A.GetZ(),B.GetX(),B.GetY(),B.GetZ()))

def roundtripFile(fname):
  '''Given a file, convert it and test for equivalence.
  This is a standalone function so it can be forked off'''
  try:
    mol = next(pybel.readfile('sdf',fname))

    #TODO TODO: move this after writing the formats to veryify addh
    #gets the same result from all formats (hint: it doesn't)
    mol.addh()

    sdftext = mol.write('sdf')
    moltext = mol.write('mol2')
    pdbtext = mol.write('pdb')
    #other formats?
    
    sdfmol = pybel.readstring('sdf',sdftext)
    if not sdfmol:
      return "Failed to convert to sdf"
    sdfmol.addh()
    
    msg = molsAreSame(mol,sdfmol)
    if msg: return 'sdfmol not equal: '+msg
      
    molmol = pybel.readstring('mol2',moltext)
    if not molmol:
      return "failed to convert to mol2"
    molmol.addh()
    
    msg = molsAreSame(mol,molmol)
    if msg: return 'molmol not equal: '+msg
    
    if True:
      pdbmol = pybel.readstring('pdb',pdbtext)
      if not pdbmol:
        return "failed to convert to pdb"      
      pdbmol.addh()    
      #dumpBoth(mol,pdbmol)
      msg = molsAreSame(mol,pdbmol)
      if msg: return 'pdbmol not equal: '+msg
      
    #roundtrip - this should have hydrogens
    fromsdf = pybel.readstring('sdf',sdfmol.write('sdf'))
    if not fromsdf:
      return 'failed to convert from sdf'
    
    msg = molsAreSame(mol,fromsdf)
    if msg: return 'fromsdf not same: '+msg

    frommol2 = pybel.readstring('sdf',molmol.write('sdf'))
    if not frommol2:
      return 'failed to convert from mol2'

    msg = molsAreSame(mol,frommol2)
    if msg: return "frommol2 not same: "+msg
       
    if True:
      frompdb = pybel.readstring('sdf',pdbmol.write('sdf'))
      if not frompdb:
        return 'failed to convert from pdb'
        
      msg = molsAreSame(mol,frompdb)
      if msg: return "frompdb not same: "+msg
            
  except Exception as e:
    traceback.print_exc()
    return str(e)
  return None


class TestSuite(unittest.TestCase):

    def setUp(self):
        self.assertTrue(ob is not None, "Failed to import the openbabel module")
        
    def canFindFile(self, filename):
        self.assertTrue(os.path.exists(filename),
                        "Cannot find the file '%s'" % filename)

    def getTestFile(self, filename):
        here = sys.path[0]
        fullpath = os.path.join(here, filename)
        self.canFindFile(fullpath)
        return fullpath    
        
    @contextmanager
    def dummy_context(self, enter_result=None):
        yield enter_result
        
    def subtest(self, enter_result=None):
        try:
            return self.subTest()
        except AttributeError: # python < 3.4
            return self.dummy_context()

            
    def testRoundtrip(self):
        """Verify PDB ligand properties are mainted through conversion"""
        root = self.getTestFile('pdb_ligands_sdf')
        #sometimes openbabel segfaults, so fork off each test        
        processes_pool = Pool(1)
        for (i,fname) in enumerate(glob.glob(os.path.join(root,'*.sdf'))):
          with self.subtest():
            try:
              result = processes_pool.apply_async(roundtripFile, args = (fname, ))
              ret = result.get(timeout=10) #surely 10 seconds is long enough
              if ret != None:
                print(i,fname)
                print(ret)
                self.assertTrue(ret == None, ret)
            except TimeoutError:
              print(i,fname)
              print("Timeout or segfault")
              self.assertTrue(False,"Timeout or segfault with %s"%fname)
          

          


if __name__ == "__main__":
  ob.obErrorLog.SetOutputLevel(ob.obError)  #ignore warnings
  if len(sys.argv) > 1:
    ret = roundtripFile(sys.argv[1])
    print(ret)
  else:
    unittest.main()
