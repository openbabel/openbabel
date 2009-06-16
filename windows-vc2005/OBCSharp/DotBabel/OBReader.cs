using System.Collections.Generic;
using System.Linq;
using System.IO;
using OpenBabel;

namespace DotBabel
{   
    public static class OBReader
    {
        private static readonly OBConversion converter;
        private static readonly VectorString obFormats;

        static OBReader()
        {
            converter = new OBConversion();
            obFormats = converter.GetSupportedInputFormat();
        }

        /// <summary>
        /// Reads a file and returns the first molecule
        /// in the file.
        /// </summary>
        /// <param name="path">path to file</param>
        /// <returns></returns>
        public static OBMol ReadMol(string path)
        {
            return ReadFiles(path).First();
        }

        /// <summary>
        /// Reads all the structures stored in all the files specified
        /// specified.
        /// </summary>
        /// <param name="fileNames">The paths to the files to be read.</param>
        /// <returns>An IEnumerable over the structures contained in the files</returns>
        public static IEnumerable<OBMol> ReadFiles(params string[] fileNames)
        {
                OBMol obMol;
                foreach (string fileName in fileNames)
                {
                    obMol = new OBMol();
                    //would be better to record the error and move on
                    //then make the class non-static and allow access to the
                    //errors
                    if (!File.Exists(fileName))
                        throw new FileNotFoundException("Could not find file : " + fileName);

                    converter.SetInFormat(OBConversion.FormatFromExt(fileName));
                    converter.ReadFile(obMol, fileName);
                    yield return obMol;

                    while (converter.Read(obMol))
                        yield return obMol;
                    
            }
        }//end ReadFiles(params string[] files)

        /// <summary>
        /// Reads all the structures stored in all the files specified
        /// specified.
        /// </summary>
        /// <param name="fileNames">The paths to the files to be read.</param>
        /// <returns>An IEnumerable over the structures contained in the files</returns>
        /// 
        public static IEnumerable<OBMol> ReadFiles(params FileInfo[] files)
        {
            OBMol obMol;

            foreach (FileInfo file in files)
            {
                obMol = new OBMol();
                //would be better to record the error and move on
                //then make the class non-static and allow access to the
                //errors
                if (!file.Exists)
                    throw new FileNotFoundException("Could not find file : " + file.Name);
                
                converter.SetInFormat(OBConversion.FormatFromExt(file.Name));
                converter.ReadFile(obMol, file.FullName);
                yield return obMol;

                while (converter.Read(obMol))
                    yield return obMol;

            }
        }

        public static VectorString Formats
        {
            get { return obFormats; }
        }

    }//end class OBReader
}
