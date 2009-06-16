using System.IO;

namespace DotBabel
{
    using OpenBabel;

    /// <summary>
    /// A facacde class for using the OBConversion object to
    /// convert files. This class can be used to convert files
    /// without exposing the rest of the OBDotNet API.
    /// </summary>
    public static class OBConvert
    {
        private static readonly OBConversion converter;
        /// <summary>
        /// static constructor to initalize the OBConversion
        /// object
        /// </summary>
        static OBConvert()
        {
            converter = new OBConversion();
        }

        public static bool ConvertFile(string inpFile, string outFile)
        {
            return ConvertFile(inpFile, outFile, false);
        }

        public static bool ConvertFile(string inpFile, string outFile, bool deleteSrcFile)
        {

            converter.SetInFormat(OBConversion.FormatFromExt(inpFile));
            converter.SetOutFormat(OBConversion.FormatFromExt(outFile));

            converter.OpenInAndOutFiles(inpFile, outFile);
            converter.Convert();
            converter.CloseOutFile();

            if (deleteSrcFile)
                File.Delete(inpFile);

            return false;
        }

        public static bool ConvertBatch(string[] fileNames, string outputType)
        {
            string outFilePath;

            foreach (string fileName in fileNames)
            {
                if (File.Exists(fileName))
                {
                    outFilePath = fileName.Substring(0, fileName.LastIndexOf(".") + 1) + outputType;
                    ConvertFile(fileName, outFilePath, false);
                }
            }

            return false;
        }


    }//end class OBConvert
}
