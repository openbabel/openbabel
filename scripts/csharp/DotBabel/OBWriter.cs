using System.Collections.Generic;
using System.Text;
using OpenBabel;

namespace DotBabel
{

    //Make sure report is in the list of allowed formats

    public static class OBWriter
    {
        private static OBConversion writer;
        private static readonly VectorString obFormats;

        static OBWriter()
        {
            writer = new OBConversion();
        }

        public static bool WriteFile(OBMol mol, string format, string fileName)
        {
            writer.SetOutFormat(format);

            return writer.WriteFile(mol, fileName);
        }

        public static bool WriteFile(OBMol mol, string fileName)
        {
            return WriteFile(mol, OBConversion.FormatFromExt(fileName).GetID(), fileName);
        }

        public static bool WriteFiles(IEnumerable<OBMol> mols, string fileName)
        {

            writer.SetOutFormat(OBConversion.FormatFromExt(fileName));

            foreach (OBMol m in mols)
            {
                writer.WriteFile((OpenBabel.OBMol)m, fileName);
            }
            return true;
        }

        public static bool WriteFile(IEnumerable<OBMol> mols, string fileName)
        {

            writer.SetOutFormat(OBConversion.FormatFromExt(fileName));

            StringBuilder outFileText = new StringBuilder(10000000);

            foreach (OBMol m in mols)
            {
                outFileText.Append(writer.WriteString((OpenBabel.OBMol)m));
            }

          
            System.IO.File.WriteAllText(fileName, outFileText.ToString());
            return true;
        }

        public static string GetFileText(OBMol mol, string format)
        {

            return writer.WriteString(mol);
        }

        public static VectorString Formats
        {
            get { return obFormats; }
        }


    } //end class OBWriter
}
