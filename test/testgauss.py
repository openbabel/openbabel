"""
Test for reading Gaussian files, determining bond orders and 
GAFF atom types. The test reads a correct SDF file with
charges and bond orders specifies, and then compares the 
molecular structure to what is generated based on the 
Gaussian output.

If the script is run with an addition (random) flag, more
verbose output will be produced.
"""
import os, sys, glob, unittest
from testbabel import BaseTest

try:
    from openbabel import openbabel as ob
except:
    print("OpenBabel not found, stopping test.")
    sys.exit(0)

debug = False

def get_mol_dict(filename, fileformat, forcefield=None):
    molecule_dict = {"molecule": {}, "atoms": {}, "bonds": {}}
    obconversion = ob.OBConversion()
    obconversion.SetInFormat(fileformat)
    obmol = ob.OBMol()

    notatend   = obconversion.ReadFile(obmol,filename)
    title      = obmol.GetTitle()
    mol_weight = obmol.GetMolWt()
    numb_atoms = obmol.NumAtoms()
    formula    = obmol.GetFormula()
    charge     = obmol.GetTotalCharge()
    
    obmol.AssignTotalChargeToAtoms(charge)
    obmol.SetAromaticPerceived(False) 

    molecule_dict["molecule"].update({"title": title, "mol_weight": mol_weight, "numb_atoms": numb_atoms, "formula": formula, "charge": charge,  "multiplicity": None })

    if forcefield:
        ff = ob.OBForceField.FindForceField(forcefield)
        if not ff:
            sys.exit("No OpenBabel support for force field %s" % forcefield)
        if not ff.Setup(obmol):
            print("Could not setup the force field %s for %s" % ( forcefield, filename))
        else:
            if not ff.GetAtomTypes(obmol):
                print("Could not get atomtypes from force field %s for %s" % ( forcefield, filename))

    # Add the atoms
    for atom in ob.OBMolAtomIter(obmol):
        index      = atom.GetIdx()
        atomtype   = "Z"
        ffatomtype = "FFAtomType"
        if forcefield and atom.HasData(ffatomtype):
            atp = atom.GetData(ffatomtype)
            if atp:
                atomtype = atp.GetValue()
        X = atom.GetX()
        Y = atom.GetY()
        Z = atom.GetZ()
        atomic_number = atom.GetAtomicNum()
        mass = atom.GetExactMass()
        molecule_dict["atoms"].update({index: {}})
        molecule_dict["atoms"][index] = {"atomic_number": atomic_number, "atomtype": atomtype, "mass": mass, "X": X, "Y": Y, "Z": Z}
    
    # Add the bonds
    for bond in ob.OBMolBondIter(obmol):
        bbb = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        molecule_dict["bonds"][bbb] = bond.GetBondOrder()
        
    return molecule_dict

def run_one(filename, forcefield, filetype):
    if not os.path.exists(filename):
        print("File %s does not exist" % ( filename ))
        return None, None, None, None
    
    moldict = get_mol_dict(filename, filetype, forcefield)
    molname = moldict["molecule"]["title"]
    atypes = []
    for ai in range(1, 1+moldict["molecule"]["numb_atoms"]):
        atype = "X"
        if ai in moldict["atoms"]:
            if moldict["atoms"][ai]["atomtype"]:
                atype = moldict["atoms"][ai]["atomtype"]
        else:
            if debug:
                print("No atom %d in %s with %d atoms" % ( ai, molname, moldict["molecule"]["numb_atoms"]))
        atypes.append(atype)
    btypes = []
    for bb in moldict["bonds"]:
        thisbond = ("%d-%d:%d" % ( bb[0], bb[1], moldict["bonds"][bb]))
        btypes.append(thisbond)
    return moldict["molecule"]["charge"], moldict["molecule"]["formula"], atypes, btypes

def atp_equal(a, b):
    # Atom type checking taking into account GAFF atom type assignmnet
    pairs = [ ( "cc", "cd" ), ( "ce", "cf" ), ( "cp", "cq" ), ( "nc", "nd" ), ( "ne", "nf" ), ( "pc", "pd" ), ( "pe", "pf" ) ]
    if a == b:
        return True
    else:
        for p in pairs:
            if (a == p[0] and b == p[1]) or (a == p[1] and b == p[0]):
                return True
    return False

class TestGauss(BaseTest):
    """Test reading files from Gaussian"""

    def compare_atom_types(self, reference, actual):
        comp = ""
        for i in range(len(reference)):
            if reference[i] != actual[i]:
                if not atp_equal(reference[i],actual[i]):
                    comp = comp + ( " %d: %s != %s" % ( i+1, reference[i], actual[i] ))
                    # We only do the assert after testing ourselves
                    # since the comparison is non-trivial
                    # When debugging we do not do the assertion, since
                    # we will not get output.
                    if not debug:
                        self.assertEqual(reference[i],actual[i])
                    
        return comp

    def compare_bond_orders(self, reference, actual):
        refhash = {}
        for r in reference.split():
            rrr = r.split(":")
            refhash[rrr[0]] = rrr[1]

        comp = ""
        for a in actual:
            aaa = a.split(":")
            if aaa[0] in refhash:
                if not debug:
                    self.assertEqual(aaa[1], refhash[aaa[0]])
                if aaa[1] != refhash[aaa[0]]:
                    comp = ( " ref %s:%s actual %s" % ( aaa[0], refhash[aaa[0]], a ))
            else:
                bbb = aaa[0].split("-")
                bnew = ( "%s-%s" % ( bbb[1], bbb[0] ))
                if bnew in refhash:
                    if not debug:
                        self.assertEqual(bnew, refhash[bnew])
                    if aaa[1] != refhash[bnew]:
                        comp += ( " ref %s:%s actual %s," % ( bnew, refhash[bnew], a ))
                else:
#                    self.assertTrue(False)
                    comp += ( " %s notfound" % bnew)
        return comp

    def compare_types(self, molname, ttype, references, actual, verbose):
        comp_atoms = ttype.find("atoms") >= 0
        if comp_atoms:
            reference = references.split()
        else:
            reference = references.split(";")[0].split()
        if len(reference) == 0:
            print("%s: no reference types for %s" % (molname, ttype))
            return
        if len(reference) != len(actual):
            extra = ""
            if verbose:
                extra = ( " ref %s actual %s" % ( reference, actual ))
            print("%s: number of %s in reference %d, actual %d%s" % ( molname, ttype, len(reference), len(actual), extra))
            return
            
        # Now we have the same, non-zero, number of atom types
        if comp_atoms:
            comp  = self.compare_atom_types(reference, actual)
        else:
            comp  = self.compare_bond_orders(references, actual)
        if len(comp) == 0:
            if debug:
                print("%s %s: Passed." % ( molname, ttype ))
            return True
        else:
            extra = ""
            if verbose:
                extra = (" ref %s actual %s" % ( reference, actual ) )
            if debug:
                print("%s %s: Failed.%s%s" % ( molname, ttype, comp, extra ) )
            return False
    
    def compare_sdf_log(self, filedir, forcefield, verbose):
        sdfs      = filedir + "/*.sdf"
        mol_list  = glob.glob(sdfs)
        filetypes = [ "sdf", "g09" ]
        summary   = { "atoms": 0, "bonds": 0 }
        passed    = True
        for mol in mol_list:
            atypes = {}
            btypes = {}
            qtot   = {}
            formula= {}
            failed = False
            for filetype in filetypes:
                if filetype == "sdf":
                    filename = mol
                else:
                    filename = mol[:-3] + "log.gz"
                qtot[filetype], formula[filetype], atypes[filetype], btypes[filetype] = run_one(filename, forcefield, filetype)
                if (qtot[filetype]    == None or 
                    formula[filetype] == None or
                    atypes[filetype]  == None or
                    btypes[filetype]  == None):
                    failed = True
                    passed = False
 
            if not failed:
                different = (qtot[filetypes[0]]    != qtot[filetypes[1]] or
                             formula[filetypes[0]] != formula[filetypes[1]])
                             
                # Now compare atom types
                ref = ""
                for at in atypes[filetypes[0]]:
                    ref += (" %s" % at)
                if not self.compare_types(mol, forcefield+"-"+"-atoms", ref, atypes[filetypes[1]], verbose):
                    summary["atoms"] += 1
                    different = True
                
                # Now compare bond types
                ref = ""
                for bt in btypes[filetypes[0]]:
                    ref += (" %s" % bt)
                if not self.compare_types(mol, forcefield+"-"+"-bonds", ref, btypes[filetypes[1]], verbose):
                    summary["bonds"] += 1
                    different = True
                # Write the atom and bond types, for sdf only
                if different:
                    passed = False
                    if debug:
                        for filetype in filetypes:
                            msg = ("%s %s qtot %s %s|" % ( mol[:-4], filetype, qtot[filetype], formula[filetype] ) )
                            for i in range(len(atypes[filetype])):
                                msg += (" %s" % atypes[filetype][i])
                            msg += ("|")
                            for i in range(len(btypes[filetype])):
                                msg += (" %s" % btypes[filetype][i])
                            print(msg)
            for ab in summary.keys():
                if debug:
                    print("%d error(s) in %s-%s" % ( summary[ab], forcefield, ab))

        return passed

    
    def testGauss(self):
        filedir = os.path.join(os.path.dirname(__file__), 'testgauss')
        result = self.compare_sdf_log(filedir, "gaff", True)
        self.assertTrue(result)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        debug = True
    testsuite = []
    for myclass in [ TestGauss ]:
        suite = unittest.TestLoader().loadTestsFromTestCase(myclass)
        testsuite.append(suite)
    unittest.TextTestRunner().run(unittest.TestSuite(testsuite))

