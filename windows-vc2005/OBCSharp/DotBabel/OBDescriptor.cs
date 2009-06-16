using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OpenBabel;

namespace DotBabel
{
    public static partial class DescCalc
    {
        private static Dictionary<string, OBDescriptor> obDesc;

        static DescCalc()
        {
            obDesc = new Dictionary<string, OBDescriptor>();
        }


        public static double EstimateProperty(OBMol mol, string propName)
        {
            if (!obDesc.Keys.Contains(propName))
                obDesc[propName] = OBDescriptor.FindType(propName);

            return obDesc[propName].Predict(mol);
        }

        public static string Estimate(OBMol mol, params string[] propNames)
        {
            return OBDescriptor.GetValues(mol, "");
        }
    }

    [Flags]
    public enum Descriptors
    { 
        ALL = 0 & 2 & 4, LOGP = 0, MR=2,TPSA=4
    }
}
