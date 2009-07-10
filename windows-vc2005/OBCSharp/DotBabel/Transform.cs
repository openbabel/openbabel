using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OpenBabel;

namespace DotBabel
{
    public delegate OBMol Transform(OBMol mol);

 
    public static class TransformOperations
    {
        public static Transform And(this Transform trans, Transform trans2)
        {
            return input => trans2(trans(input));
        }

        public static Transform Or(this Transform trans, Transform trans2)
        {
            return input => trans(input) ?? trans2(input);
        }
    }

    public class BasicTransforms
    {
        public static Transform SmartsReplace(string smarts, string repSmiles)
        {
            return Replace(a => a.MatchesSMARTS(smarts), repSmiles);
            
        }//end SmartsReplace


        public static Transform Replace(Func<OBAtom,bool> selector, string smi)
        {
           // uint n = input.Atoms().First(selector).GetIdx();
            return null;
        }

        public static Transform Replace(int idx, string smi)
        {
            return Replace(idx, Smiles.ParseSmiles(smi));
        }

        public static Transform Replace(int idx, OBMol mol)
        {
            return input =>
            {
                OBMol retMol = new OBMol(input);

                OBAtom target = input.GetAtom(idx);
                uint joinIdx = 0;

             
                foreach(OBAtom a in mol.Atoms())
                {
                    mol.AddAtom(a);
                }

                foreach (OBAtom neighbor in target.Neighbors())
                { 
                   
                }

                return retMol;
            };
        }
    }

  
}