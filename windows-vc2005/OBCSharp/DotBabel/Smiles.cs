using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OpenBabel;

namespace DotBabel
{
    public static class Smiles
    {
        private static readonly OBConversion smilesReader;

        static Smiles()
        {
            smilesReader = new OBConversion();
            smilesReader.SetInAndOutFormats("smi", "can");
        }

        public static OBMol ParseSmiles(string smiles)
        {   
            OBMol mol = new OBMol();
            smilesReader.ReadString(mol, smiles);
            return mol;
        }

        public static string GenSmiles(OBMol mol)
        {
            
            return smilesReader.WriteString(mol,true);
        }
    }
}
