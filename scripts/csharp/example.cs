using System;
// using System.Collections.Generic;
// using System.Text;
using OpenBabel;

namespace MyConsoleApplication
{
    class Program
    {
        static void Main(string[] args)
        {
            OBConversion obconv = new OBConversion();
            obconv.SetInFormat("smi");
            OBMol mol = new OBMol();
            obconv.ReadString(mol, "CCC");
            System.Console.WriteLine(mol.GetMolWt());
            //System.Console.ReadKey();
        }
    }
}

