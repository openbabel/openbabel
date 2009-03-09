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
            return OBWriter.GetFileText(mol, "can");
        }

        public static string AsInChI(this OBMol mol)
        {
            return OBWriter.GetFileText(mol, "inchi");
        }
    }
}
