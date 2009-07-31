using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OpenBabel;

namespace DotBabel
{
    public static class OBMolExtensions
    {
        public static double CalcDescriptor(this OBMol mol, string descType)
        {
            return DescCalc.EstimateProperty(mol, descType);
        }

        public static string AsSmiles(this OBMol mol)
        {
            return Smiles.GenSmiles(mol);
        }

        public static string AsInChI(this OBMol mol)
        {
            return OBWriter.GetFileText(mol, "inchi");
        }

        public static bool Make3D(this OBMol mol)
        {
            OBOp Gen3D = OBOp.FindType("Gen3D");
            return Gen3D.Do(mol);
        }

        public static void AddAtoms(this OBMol mol,IEnumerable<OBAtom> atoms)
        {
            foreach (OBAtom a in atoms)
                mol.AddAtom(a);
        }

        public static void AddBonds(this OBMol mol, IEnumerable<OBBond> bonds)
        {
            foreach (OBBond b in bonds)
                mol.AddBond(b);
        }

        public static OBMol DoTransformation(this OBMol mol,Transform tr)
        {
            return tr(mol);
        }
    }
}
